"""Classes representing a workflow."""

import sys
import os
import os.path
import string
import subprocess

import gwf_workflow
from gwf_workflow.jobs import JOBS_QUEUE


def _escape_file_name(filename):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in filename if c in valid_chars)


def _escape_job_name(job_name):
    valid_chars = "_%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in job_name if c in valid_chars)


_remembered_files = {}  # cache to avoid too much stat'ing
_remembered_timestamps = {}  # Use this to avoid too many stats that slows down the script


def _file_exists(filename):
    if filename not in _remembered_files:
        _remembered_files[filename] = os.path.exists(filename)
    return _remembered_files[filename]


def _get_file_timestamp(filename):
    if filename not in _remembered_timestamps:
        _remembered_timestamps[filename] = os.path.getmtime(filename)
    return _remembered_timestamps[filename]


def _make_absolute_path(working_dir, filename):
    if os.path.isabs(filename):
        abspath = filename
    else:
        abspath = os.path.join(working_dir, filename)
    return os.path.normpath(abspath)


class Node(object):
    """Class handling targets. Stores the info for executing them."""

    def __init__(self, target):
        self.target = target
        self.depends_on = set()
        self.dependents = set()

        self.script_dir = _make_absolute_path(self.target.working_dir, '.scripts')
        self.script_name = _make_absolute_path(self.script_dir, _escape_file_name(self.target.name))

        self.cached_should_run = None
        self.reason_to_run = None

    @property
    def should_run(self):
        """Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. """

        if self.cached_should_run is not None:
            return self.cached_should_run

        if len(self.target.options['output']) == 0:
            self.reason_to_run = 'Sinks (targets without output) should always run'
            self.cached_should_run = True
            return True

        for outfile in self.target.options['output']:
            if not _file_exists(_make_absolute_path(self.target.working_dir, outfile)):
                self.reason_to_run = 'Output file %s is missing' % outfile
                self.cached_should_run = True
                return True

        for infile in self.target.options['input']:
            if not _file_exists(_make_absolute_path(self.target.working_dir, infile)):
                self.reason_to_run = 'Input file %s is missing' % infile
                self.cached_should_run = True
                return True

        # If no file is missing, it comes down to the time stamps. If we
        # only have output and no input, we assume the output is up to
        # date. Touching files and adding input can fix this behaviour
        # from the user side but if we have a program that just creates
        # files we don't want to run it whenever someone needs that
        # output just because we don't have time stamped input.

        if len(self.target.options['input']) == 0:
            self.reason_to_run = "We shouldn't run"
            self.cached_should_run = False
            return False

        # if we have both input and output files, check time stamps

        youngest_in_timestamp = None
        youngest_in_filename = None
        for infile in self.target.options['input']:
            timestamp = _get_file_timestamp(_make_absolute_path(self.target.working_dir, infile))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_filename = infile
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        oldest_out_timestamp = None
        oldest_out_filename = None
        for outfile in self.target.options['output']:
            timestamp = _get_file_timestamp(_make_absolute_path(self.target.working_dir, outfile))
            if oldest_out_timestamp is None \
                    or oldest_out_timestamp > timestamp:
                oldest_out_filename = outfile
                oldest_out_timestamp = timestamp
        assert oldest_out_timestamp is not None

        # The youngest in should be older than the oldest out
        if youngest_in_timestamp >= oldest_out_timestamp:
            # we have a younger in file than an outfile
            self.reason_to_run = 'Infile %s is younger than outfile %s' % (youngest_in_filename, oldest_out_filename)
            self.cached_should_run = True
            return True
        else:
            self.reason_to_run = 'Youngest infile %s is older than ' \
                                 'the oldest outfile %s' % \
                                 (youngest_in_filename, oldest_out_filename)
            self.cached_should_run = False
            return False

    def make_script_dir(self):
        script_dir = self.script_dir
        # Don't use the _file_exists() function here. It caches its status and that won't work for the script dir.
        if not os.path.exists(script_dir):
            os.makedirs(script_dir)

    def write_script(self):
        """Write the code to a script that can be executed."""

        self.make_script_dir()
        f = open(self.script_name, 'w')

        print >> f, "#!/bin/bash"

        from gwf_workflow import BACKEND
        BACKEND.write_script_header(f, self.target.options)
        print >> f,

        print >> f, '# GWF generated code ...'
        print >> f, 'cd %s' % self.target.working_dir
        BACKEND.write_script_variables(f, self.target.options)
        print >> f

        print >> f, '# Script from workflow'
        print >> f, self.target.spec

    @property
    def job_in_queue(self):
        return JOBS_QUEUE.get_database(self.target.working_dir).in_queue(self.target.name)

    @property
    def job_queue_status(self):
        return JOBS_QUEUE.get_database(self.target.working_dir).get_job_status(self.target.name)

    @property
    def job_id(self):
        if self.job_in_queue:
            return JOBS_QUEUE.get_database(self.target.working_dir).get_job_id(self.target.name)
        else:
            return None

    def set_job_id(self, job_id):
        JOBS_QUEUE.get_database(self.target.working_dir).set_job_id(self.target.name, job_id)

    def submit(self, dependents):
        if self.job_in_queue:
            return self.job_id

        self.write_script()
        command = ['qsub', '-N', self.target.name]
        dependents_ids = [dependent.job_id for dependent in dependents]
        if len(dependents_ids) > 0:
            command.append('-W')
            command.append('depend=afterok:{}'.format(':'.join(dependents_ids)))
        command.append(self.script_name)

        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        self.set_job_id(job_id)
        return job_id

    def get_existing_outfiles(self):
        """Get list of output files that already exists."""
        result = []
        for outfile in self.target.options['output']:
            filename = _make_absolute_path(self.target.working_dir, outfile)
            if _file_exists(filename):
                result.append(filename)
        return result

    def clean_target(self):
        """Delete all existing outfiles."""
        for fname in self.get_existing_outfiles():
            os.remove(fname)

    def __str__(self):
        return str(self.target)

    __repr__ = __str__  # not really the correct use of __repr__ but easy
    # for printing output when testing...


def schedule(nodes, target_name):
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
            scheduled.add(node.target.name)

        processed.add(node)

    dfs(root)

    return job_schedule, scheduled


# # WRAPPING IT ALL UP IN A WORKFLOW...
class Workflow:
    """Class representing a workflow."""

    def __init__(self, targets):
        self.targets = targets

    def get_execution_schedule(self, target_name):
        execution_schedule, scheduled_tasks = schedule(self.targets, target_name)
        return execution_schedule, scheduled_tasks


def build_workflow():
    """Collect all the targets and build up their dependencies."""

    nodes = {}
    providing = {}
    for target in gwf_workflow.ALL_TARGETS.values():
        target.options['input'] = [_make_absolute_path(target.working_dir, infile)
                                   for infile in target.options['input']]
        target.options['output'] = [_make_absolute_path(target.working_dir, outfile)
                                    for outfile in target.options['output']]
        node = Node(target)
        for outfile in target.options['output']:
            if outfile in providing:
                print 'Warning: File', outfile, 'is provided by both',
                print target.name, 'and', providing[outfile].target.name
            providing[outfile] = node
        nodes[target.name] = node

    for node in nodes.values():
        for infile in node.target.options['input']:
            if infile in providing:
                provider = providing[infile]
                node.depends_on.add(provider)
                provider.dependents.add(node)
            else:
                if not _file_exists(infile):
                    print 'Target', node.target.name, 'needs file',
                    print infile, 'which does not exists and is not constructed.'
                    print 'Aborting'
                    sys.exit(2)

    return Workflow(nodes)
