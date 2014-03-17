"""Classes representing a workflow."""

import sys
import os
import os.path
import re
import string
import subprocess

import gwf_workflow


def _escape_file_name(filename):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in filename if c in valid_chars)


def _escape_job_name(job_name):
    valid_chars = "_%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in job_name if c in valid_chars)


def _file_exists(filename):
    return os.path.exists(filename)


def _get_file_timestamp(filename):
    return os.path.getmtime(filename)


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

    @property
    def should_run(self):
        """Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. """

        if len(self.target.output) == 0:
            self.reason_to_run = \
                'Sinks (targets without output) should always run'
            return True  # If we don't provide output, assume we always
            # need to run.

        for outfile in self.target.output:
            if not _file_exists(_make_absolute_path(self.target.working_dir, outfile)):
                self.reason_to_run = \
                    'Output file %s is missing' % outfile
                return True

        for infile in self.target.input:
            if not _file_exists(_make_absolute_path(self.target.working_dir, infile)):
                self.reason_to_run = \
                    'Input file %s is missing' % infile
                return True

        # If no file is missing, it comes down to the time stamps. If we
        # only have output and no input, we assume the output is up to
        # date. Touching files and adding input can fix this behaviour
        # from the user side but if we have a program that just creates
        # files we don't want to run it whenever someone needs that
        # output just because we don't have time stamped input.

        if len(self.target.input) == 0:
            self.reason_to_run = "We shouldn't run"
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
            self.reason_to_run = 'Infile %s is younger than outfile %s' % \
                                 (youngest_in_filename, oldest_out_filename)
            return True
        else:
            self.reason_to_run = 'Youngest infile %s is older than ' \
                                 'the oldest outfile %s' % \
                                 (youngest_in_filename, oldest_out_filename)
            return False

        assert False, "We shouldn't get here"


    @property
    def script_dir(self):
        return _make_absolute_path(self.target.working_dir, '.scripts')

    @property
    def jobs_dir(self):
        return _make_absolute_path(self.target.working_dir, '.jobs')

    @property
    def script_name(self):
        # Escape name to make a file name...
        escaped_name = _escape_file_name(self.target.name)
        return _make_absolute_path(self.script_dir, escaped_name)

    @property
    def job_name(self):
        # Escape name to make a file name...
        escaped_name = _escape_file_name(self.target.name)
        return _make_absolute_path(self.jobs_dir, escaped_name)

    def make_script_dir(self):
        script_dir = self.script_dir
        if _file_exists(script_dir):
            return
        os.makedirs(script_dir)

    def make_jobs_dir(self):
        jobs_dir = self.jobs_dir
        if _file_exists(jobs_dir):
            return
        os.makedirs(jobs_dir)

    def write_script(self):
        '''Write the code to a script that can be executed.'''

        self.make_script_dir()
        self.make_jobs_dir()

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
        if not _file_exists(self.job_name):
            return False
        else:
            stat = subprocess.Popen(['qstat', '-f', self.jobID],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
            for line in stat.stdout:
                line = line.strip()
                if line.startswith('job_state'):
                    self.JOB_QUEUE_STATUS = line.split()[2]
                    if self.JOB_QUEUE_STATUS == 'E':
                        # We don't consider a failed job as being
                        # in the queue
                        return False
                    else:
                        return True
            return False

    @property
    def job_queue_status(self):
        '''Get the job status if the job is in the queue.'''
        # First check if it is cached
        if hasattr(self, 'JOB_QUEUE_STATUS'):
            return self.JOB_QUEUE_STATUS

        # If it isn't, get the job status implicitly by checking
        # its queue status
        self.job_in_queue
        # and now return if it that worked, or just return None
        if hasattr(self, 'JOB_QUEUE_STATUS'):
            return self.JOB_QUEUE_STATUS
        else:
            return None

    @property
    def jobID(self):
        if _file_exists(self.job_name):
            return open(self.job_name).read().strip()
        else:
            return None

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

    print [n.target.name for n in nodes.values() if len(n.depends_on) == 0]
    print [n.target.name for n in nodes.values() if len(n.dependents) == 0]

## WRAPPING IT ALL UP IN A WORKFLOW...
class Workflow:
    '''Class representing a workflow.'''

    def __init__(self, lists, templates, targets, template_targets, wd):
        self.lists = lists
        self.templates = templates
        self.targets = targets
        self.template_targets = template_targets
        self.working_dir = wd

        # handle list transformation...
        for cmd in self.lists.values():
            if isinstance(cmd, Transform):
                input_list_name = cmd.input_list
                if input_list_name not in self.lists:
                    print "Transformation list %s uses input list %s the doesn't exist." % \
                          (cmd.name, input_list_name)
                    sys.exit(2)

                input_list = self.lists[input_list_name]

                cmd.elements = [re.sub(cmd.match_pattern, cmd.subs_pattern, input)
                                for input in input_list.elements]

        # handle list expansions in other lists
        for cmd in self.lists.values():
            def expand_lists(lst):
                new_list = []
                for elm in lst:
                    if elm.startswith('@'):
                        listname = elm[1:]
                        if listname not in self.lists:
                            print 'List %s references unknown list %s.' % \
                                  (cmd.name, listname)
                            sys.exit(2)
                        new_list.extend(self.lists[listname].elements)
                    else:
                        new_list.append(elm)
                return new_list

            cmd.elements = expand_lists(cmd.elements)


        # handle templates and template instantiations
        for name, tt in self.template_targets.items():
            for target_code in tt.instantiate_target_code(self):
                target = parser.parse_target(target_code, tt.working_dir)
                if target.name in self.targets:
                    print 'Instantiated template %s has the same name as an existing target' % \
                          target.name
                    sys.exit(2)
                self.targets[target.name] = target

        # expand lists in input and output lists for the targets
        for target in self.targets.values():
            def expand_lists(lst):
                new_list = []
                for elm in lst:
                    if elm.startswith('@'):
                        listname = elm[1:]
                        if listname not in self.lists:
                            print 'Target %s references unknown list %s.' % \
                                  (target.name, listname)
                            sys.exit(2)
                        new_list.extend(self.lists[listname].elements)
                    else:
                        new_list.append(elm)
                return new_list

            target.input = expand_lists(target.input)
            target.output = expand_lists(target.output)

            # make all files absolute and normalised so different ways of
            # referring to the same file actually works.
            # For obvious reasons this has to go after list expansion...
            target.input = [_make_absolute_path(target.working_dir, fname)
                            for fname in target.input]
            target.output = [_make_absolute_path(target.working_dir, fname)
                             for fname in target.output]


        # collect the output files so we know who can build them.
        self.providers = dict()
        for target in self.targets.values():
            for output_file in target.output:
                assert output_file not in self.providers
                self.providers[output_file] = target

        # now get dependencies for each target...
        for target in self.targets.values():
            dependencies = []
            for input_file in target.input:
                if input_file in self.providers:
                    dependencies.append((input_file,
                                         self.providers[input_file]))
                else:
                    sysfile = SystemFile(input_file, self.working_dir)
                    dependencies.append((input_file, sysfile))
            target.dependencies = dependencies

        # build the dependency graph    
        self.dependency_graph = DependencyGraph(self)


    def get_submission_script(self, target_name):
        '''Generate the script used to submit the tasks.'''

        target = self.targets[target_name]
        schedule, scheduled_tasks = self.dependency_graph.schedule(target.name)

        script_commands = []
        for job in schedule:

            # skip dummy tasks that we shouldn't submit...
            if job.task.dummy:
                continue

            if not job.task.can_execute:
                print job.task.execution_error
                import sys;

                sys.exit(2)

            # If the job is already in the queue, just get the ID
            # into the shell command used later for dependencies...
            if job.task.job_in_queue:
                command = ' '.join([
                    '%s=`' % job.name,
                    'cat', job.task.job_name,
                    '`'])
                script_commands.append(command)

            else:
                # make sure we have the scripts for the jobs we want to
                # execute!
                job.task.write_script()

                dependent_tasks = set(node.name
                                      for _, node in job.dependencies
                                      if node.name in scheduled_tasks)
                if len(dependent_tasks) > 0:
                    depend = '-W depend=afterok:$%s' % \
                             ':$'.join(dependent_tasks)
                else:
                    depend = ''

                script = job.task.script_name
                command = ' '.join([
                    '%s=`' % job.name,
                    'qsub -N %s' % job.name,
                    depend,
                    script,
                    '`'])
                script_commands.append(command)
                script_commands.append(' '.join([
                    'echo', ('$%s' % job.name), '>', job.task.job_name]))

        return '\n'.join(script_commands)


    def get_local_execution_script(self, target_name):
        '''Generate the script needed to execute a target locally.'''

        target = self.targets[target_name]
        schedule, scheduled_tasks = self.dependency_graph.schedule(target.name)

        script_commands = []
        for job in schedule:

            # skip dummy tasks that we shouldn't submit...
            if job.task.dummy:
                continue

            if not job.task.can_execute:
                print job.task.execution_error
                import sys;

                sys.exit(2)

            # make sure we have the scripts for the jobs we want to
            # execute!
            job.task.write_script()

            script = open(job.task.script_name, 'r').read()
            script_commands.append('# computing %s' % job.name)
            script_commands.append(script)

        return '\n'.join(script_commands)


