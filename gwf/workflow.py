from __future__ import absolute_import
from __future__ import print_function

import os.path
import sys

import gwf

from gwf.colours import *
from gwf.helpers import make_list, make_absolute_path, file_exists, get_file_timestamp
from gwf.jobs import JOBS_QUEUE


# This global variable will hold all the targets after a workflow script has completed.
# gwf will use this list for its further processing.
ALL_TARGETS = {}


class Target(object):
    """Class handling targets. Stores the info for executing them."""

    def __init__(self, name, options, spec):
        self.name = name
        self.spec = spec

        self.options = {
            'input': [],
            'output': [],
            'nodes': 1,
            'cores': 1,
            'memory': "4g",
            'walltime': '120:00:00',
        }

        known_options_without_defaults = {'queue', 'account', 'constraint', 'mail_user', 'mail_type'}

        for k in options.keys():
            if k in self.options or k in known_options_without_defaults:
                self.options[k] = options[k]
            else:
                print('Warning:, Target', self.name, 'has unknown option', k)

        # handle that input and output can be both lists and single file names
        self.options['input'] = filter(None, make_list(self.options['input']))
        self.options['output'] = filter(None, make_list(self.options['output']))

        self.depends_on = set()
        self.dependents = set()

        self.cached_node_should_run = None
        self.reason_to_run = None
        self.cached_should_run = None

    @property
    def node_should_run(self):
        """Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. """

        if self.cached_node_should_run is not None:
            return self.cached_node_should_run

        if len(self.options['output']) == 0:
            self.reason_to_run = 'Sinks (targets without output) should always run'
            self.cached_node_should_run = True
            return True

        for outfile in self.options['output']:
            if not file_exists(make_absolute_path(gwf.WORKING_DIR, outfile)):
                self.reason_to_run = 'Output file %s is missing' % outfile
                self.cached_node_should_run = True
                return True

        for infile in self.options['input']:
            if not file_exists(make_absolute_path(gwf.WORKING_DIR, infile)):
                self.reason_to_run = 'Input file %s is missing' % infile
                self.cached_node_should_run = True
                return True

        # If no file is missing, it comes down to the time stamps. If we
        # only have output and no input, we assume the output is up to
        # date. Touching files and adding input can fix this behaviour
        # from the user side but if we have a program that just creates
        # files we don't want to run it whenever someone needs that
        # output just because we don't have time stamped input.

        if len(self.options['input']) == 0:
            self.reason_to_run = "We shouldn't run"
            self.cached_node_should_run = False
            return False

        # if we have both input and output files, check time stamps

        youngest_in_timestamp = None
        youngest_in_filename = None
        for infile in self.options['input']:
            timestamp = get_file_timestamp(make_absolute_path(gwf.WORKING_DIR, infile))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_filename = infile
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        oldest_out_timestamp = None
        oldest_out_filename = None
        for outfile in self.options['output']:
            timestamp = get_file_timestamp(make_absolute_path(gwf.WORKING_DIR, outfile))
            if oldest_out_timestamp is None \
                    or oldest_out_timestamp > timestamp:
                oldest_out_filename = outfile
                oldest_out_timestamp = timestamp
        assert oldest_out_timestamp is not None

        # The youngest in should be older than the oldest out
        if youngest_in_timestamp > oldest_out_timestamp:
            # we have a younger in file than an outfile
            self.reason_to_run = 'Infile %s is younger than outfile %s' % (youngest_in_filename, oldest_out_filename)
            self.cached_node_should_run = True
            return True
        else:
            self.reason_to_run = 'Youngest infile %s is older than ' \
                                 'the oldest outfile %s' % \
                                 (youngest_in_filename, oldest_out_filename)
            self.cached_node_should_run = False
            return False

    @property
    def should_run(self):
        if self.cached_should_run is not None:
            return self.cached_should_run

        # If this target should run, or any of the targets and this target depends on should run, this node should run.
        self.cached_should_run = self.node_should_run or any(n.should_run for n in self.depends_on)
        return self.cached_should_run

    @property
    def job_in_queue(self):
        return JOBS_QUEUE.get_database(gwf.WORKING_DIR).in_queue(self.name)

    @property
    def job_queue_status(self):
        return JOBS_QUEUE.get_database(gwf.WORKING_DIR).get_job_status(self.name)

    @property
    def job_id(self):
        return JOBS_QUEUE.get_database(gwf.WORkING_DIR).get_job_id(self.name)

    def set_job_id(self, job_id):
        JOBS_QUEUE.get_database(gwf.WORKING_DIR).set_job_id(self.name, job_id)

    def submit(self, dependents):
        if self.job_in_queue:
            return self.job_id

        dependents_ids = [dependent.job_id for dependent in dependents]
        try:
            job_id = gwf.backends.BACKEND.submit_command(self, self.script_name, dependents_ids)
            self.set_job_id(job_id)
        except OSError as ex:
            print()
            print(COLORS['red'], COLORS['bold'])
            print('ERROR:', CLEAR, end=' ')
            print("Couldn't execute the submission command {}'{}'{}.".format(COLORS['bold'], ' '.join(self.script_name), CLEAR))
            print(ex)
            print(COLORS['red'])
            print("Quiting submissions", CLEAR)
            print()

            sys.exit(2)

        return job_id

    def get_existing_outfiles(self):
        """Get list of output files that already exists."""
        result = []
        for outfile in self.options['output']:
            filename = make_absolute_path(gwf.WORKING_DIR, outfile)
            if file_exists(filename):
                result.append(filename)
        return result

    def clean_target(self):
        """Delete all existing outfiles."""
        for fname in self.get_existing_outfiles():
            os.remove(fname)

    def __str__(self):
        return str(self)

    __repr__ = __str__  # not really the correct use of __repr__ but easy
    # for printing output when testing...


def dependencies(nodes, target_name):
    """Return all tasks necessary for building the target.

    The set of tasks is just returned as set.
    """
    root = nodes[target_name]

    # Working with a list to preserve the order. It makes lookups slower but hopefully these sets
    # won't be terribly long ... if it becomes a problem it is easy enough to fix it.
    processed = []

    def dfs(node):
        if node in processed:
            return
        else:
            for dep in node.depends_on:
                dfs(dep)
            processed.append(node)

    dfs(root)
    return processed


def get_execution_schedule(nodes, target_name):
    """Linearize the targets to be run.

    Returns a list of tasks to be run (in the order they should run or
    be submitted to the cluster to make sure dependencies are handled
    correctly) and a set of the names of tasks that will be scheduled
    (to make sure dependency flags are set in the submission command).
    """

    root = nodes[target_name]

    # If the target is already in the queue we just dismiss the scheduling
    # right away... this because we need to handle dependent nodes in the
    # queue differently, since for those we need wait for completion.
    if root.job_in_queue:
        return [], set()

    processed = set()
    scheduled = set()
    job_schedule = []

    def dfs(node):
        if node in processed:
            # we have already processed the node, and
            # if we should run the target name is scheduled
            # otherwise it isn't.
            return node.target.name in scheduled

        # schedule all dependencies before we schedule this task
        for dep in node.depends_on:
            dfs(dep)

        # If this task needs to run, then schedule it
        if node.job_in_queue or node.should_run:
            job_schedule.append(node)
            scheduled.add(node.name)

        processed.add(node)

    dfs(root)

    return job_schedule, scheduled


def build_workflow():
    """Collect all the targets and build up their dependencies."""

    nodes = {}
    providing = {}
    for target in ALL_TARGETS.values():
        target.options['input'] = [make_absolute_path(gwf.WORKING_DIR, infile)
                                   for infile in target.options['input']]
        target.options['output'] = [make_absolute_path(gwf.WORKING_DIR, outfile)
                                    for outfile in target.options['output']]
        #node = Node(target)
        for outfile in target.options['output']:
            if outfile in providing:
                print('Warning: File', outfile, 'is provided by both', end=' ')
                print(target.name, 'and', providing[outfile].target.name)
            providing[outfile] = target
        nodes[target.name] = target

    for node in nodes.values():
        for infile in node.options['input']:
            if infile in providing:
                provider = providing[infile]
                node.depends_on.add(provider)
                provider.dependents.add(node)
            else:
                if not file_exists(infile):
                    print('Target', node.target.name, 'needs file', end=' ')
                    print(infile, 'which does not exists and is not constructed.')
                    print('Aborting')
                    sys.exit(2)

    return nodes


def split_tasks(tasks):
    up_to_date, in_queue, to_schedule = [], [], []
    for task in tasks:
        if task.job_in_queue:
            in_queue.append(task)
        elif task.should_run:
            to_schedule.append(task)
        else:
            up_to_date.append(task)
    return up_to_date, in_queue, to_schedule