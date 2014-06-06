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
def _file_exists(filename):
    if filename not in _remembered_files:
        _remembered_files[filename] = os.path.exists(filename)
    return _remembered_files[filename]

_remembered_timestamps = {}  # Use this to avoid too many stats that slows down the script
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

        if len(self.target.output) == 0:
            self.reason_to_run = 'Sinks (targets without output) should always run'
            self.cached_should_run = True
            return True

        for outfile in self.target.output:
            if not _file_exists(_make_absolute_path(self.target.working_dir, outfile)):
                self.reason_to_run = 'Output file %s is missing' % outfile
                self.cached_should_run = True
                return True

        for infile in self.target.input:
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

        if len(self.target.input) == 0:
            self.reason_to_run = "We shouldn't run"
            self.cached_should_run = False
            return False

        # if we have both input and output files, check time stamps

        youngest_in_timestamp = None
        youngest_in_filename = None
        for infile in self.target.input:
            timestamp = _get_file_timestamp(_make_absolute_path(self.target.working_dir, infile))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_filename = infile
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        oldest_out_timestamp = None
        oldest_out_filename = None
        for outfile in self.target.output:
            timestamp = _get_file_timestamp(_make_absolute_path(self.target.working_dir, outfile))
            if oldest_out_timestamp is None \
                    or oldest_out_timestamp > timestamp:
                oldest_out_filename = outfile
                oldest_out_timestamp = timestamp
        assert oldest_out_timestamp is not None

        # The youngest in should be older than the oldest out
        if youngest_in_timestamp >= oldest_out_timestamp:
            # we have a younger in file than an outfile
            self.reason_to_run = 'Infile %s is younger than outfile %s' %  (youngest_in_filename, oldest_out_filename)
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
        if _file_exists(script_dir):
            return
        os.makedirs(script_dir)

    def write_script(self):
        """Write the code to a script that can be executed."""

        self.make_script_dir()
        f = open(self.script_name, 'w')

        # Put PBS options at the top
        for options in self.target.pbs:
            print >> f, '#PBS', options
        print >> f

        print >> f, '# GWF generated code ...'
        print >> f, 'cd %s' % self.target.working_dir
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
        command = ['qsub', '-N', self.target.name, self.script_name]
        dependents_ids = [dependent.job_id for dependent in dependents]
        if len(dependents_ids) > 0:
            command.append('-W depend=afterok:{}'.format(':'.join(dependents_ids)))

        qsub = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        job_id = qsub.stdout.read().strip()
        self.set_job_id(job_id)
        return job_id


    def get_existing_outfiles(self):
        """Get list of output files that already exists."""
        result = []
        for outfile in self.target.output:
            filename = _make_absolute_path(self.target.working_dir, outfile)
            if _file_exists(filename):
                result.append(filename)
        return result

    def clean_target(self):
        '''Delete all existing outfiles.'''
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
        (to make sure dependency flags are set in the qsub command).

        """

        root = nodes[target_name]

        processed = set()
        scheduled = set()
        schedule = []

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
                schedule.append(node)
                scheduled.add(node.target.name)

            processed.add(node)

        dfs(root)

        return schedule, scheduled


## WRAPPING IT ALL UP IN A WORKFLOW...
class Workflow:
    """Class representing a workflow."""

    def __init__(self, targets):
        self.targets = targets

    def get_execution_schedule(self, target_name):
        execution_schedule, scheduled_tasks = schedule(self.targets, target_name)
        return execution_schedule, scheduled_tasks

    def get_submission_script(self, target_name):
        """Generate the script used to submit the tasks."""

        execution_schedule, scheduled_tasks = schedule(self.targets, target_name)

        script_commands = []
        for job in execution_schedule:

            # If the job is already in the queue, just get the ID
            # into the shell command used later for dependencies...
            if job.job_in_queue:
                command = ' '.join([
                    '%s=`' % job.target.name,
                    'cat', job.job_name,
                    '`'])
                script_commands.append(command)

            else:
                # make sure we have the scripts for the jobs we want to
                # execute!
                job.write_script()

                dependent_tasks = set(node.target.name
                                      for node in job.depends_on
                                      if node.target.name in scheduled_tasks)
                if len(dependent_tasks) > 0:
                    depend = '-W depend=afterok:$%s' % \
                             ':$'.join(dependent_tasks)
                else:
                    depend = ''

                script = job.script_name
                command = ' '.join([
                    '%s=`' % job.target.name,
                    'qsub -N %s' % job.target.name,
                    depend,
                    script,
                    '`'])
                script_commands.append(command)
                script_commands.append(' '.join([
                    'echo', ('$%s' % job.target.name), '>', job.job_name]))

        return '\n'.join(script_commands)


    def get_local_execution_script(self, target_name):
        '''Generate the script needed to execute a target locally.'''

        execution_schedule, scheduled_tasks = schedule(self.targets, target_name)

        script_commands = []
        for job in execution_schedule:

            # make sure we have the scripts for the jobs we want to
            # execute!
            job.write_script()

            script = open(job.script_name, 'r').read()
            script_commands.append('# computing %s' % job.target.name)
            script_commands.append(script)

        return '\n'.join(script_commands)


def build_workflow():
    """Collect all the targets and build up their dependencies."""

    nodes = {}
    providing = {}
    for target in gwf_workflow.ALL_TARGETS.values():
        target.input = [_make_absolute_path(target.working_dir, infile)
                        for infile in target.input]
        target.output = [_make_absolute_path(target.working_dir, outfile)
                         for outfile in target.output]
        node = Node(target)
        for outfile in target.output:
            if outfile in providing:
                print 'Warning: File', outfile, 'is provided by both',
                print target.name, 'and', providing[outfile].target.name
            providing[outfile] = node
        nodes[target.name] = node

    for node in nodes.values():
        for infile in node.target.input:
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
