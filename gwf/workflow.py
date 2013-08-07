'''Classes representing a workflow.'''

import sys
import os, os.path
import time
import re
import string
import subprocess
from exceptions import NotImplementedError
from dependency_graph import DependencyGraph
import parser # need this to re-parse instantiated templates

def _escape_file_name(fname):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in fname if c in valid_chars)
    
def _escape_job_name(jobname):
    valid_chars = "_%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in jobname if c in valid_chars)

def _file_exists(fname):
    return os.path.exists(fname)

def _get_file_timestamp(fname):
    return os.path.getmtime(fname)
    
def _make_absolute_path(working_dir, fname):
    if os.path.isabs(fname):
        abspath = fname
    else:
        abspath = os.path.join(working_dir, fname)
    return os.path.normpath(abspath)

## TEMPLATES 
class Template:
    def __init__(self, name, wd, parameters, template):
        self.name = name
        self.working_dir = wd
        self.parameters = parameters
        self.template = template
        
    def __str__(self):
        return '@template %s %s [%s...]' % (
            self.name,
            self.parameters,
            self.template[:80].replace('\n',' ')
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...
        
    
class TemplateTarget:
    def __init__(self, name, wd, template, parameter_assignments):
        self.name = name
        self.working_dir = wd
        self.template = template
        self.assignments = parameter_assignments
    
    def instantiate_target_code(self, workflow):
        '''Instantiate a target from the template. Uses "workflow" to access
        variable names like lists and templates.'''

        if self.template not in workflow.templates:
            print 'Template-target %s refers to unknown template %s' % \
                (self.name, self.template)
            sys.exit(2)        
        template_code = workflow.templates[self.template].template
        
        def instance(assignments):
            return 'target %s\n%s' % (_escape_job_name(self.name.format(**assignments)),
                                      template_code.format(**assignments))

        # If there are variables in the instantiation we must expand them 
        # here. If we are dealing with lists we get more than one target.

        has_variables = any(v.startswith('@')
                            for v in self.assignments.values())
                            
        if has_variables:
            # We have variables, so first we try to expand them.
            expanded = dict()
            for key,val in self.assignments.items():
                if val.startswith('@'):
                    listname = val[1:]
                    if listname not in workflow.lists:
                        print 'Template target %s refers unknown list %s' %\
                            (self.name, listname)
                        sys.exit(2)
                    elements = workflow.lists[listname].elements
                    expanded[val] = elements
            
            # Make sure the expanded lists have the same length...
            lengths = set(len(elms) for elms in expanded.values())
            if len(lengths) != 1:
                print 'The lists used in target template', self.name,
                print 'have different length.'
                for name,elms in expanded.items():
                    print name, '=', elms
                sys.exit(2)
            n = lengths.pop()
            
            new_assignments = dict()
            for key,val in self.assignments.items():
                if val in expanded:
                    new_assignments[key] = expanded[val]
                else:
                    new_assignments[key] = [val] * n
            
            targets = []
            for i in xrange(n):
                inst_assignments = \
                    dict((k,v[i]) for k,v in new_assignments.items())
                targets.append(instance(inst_assignments))
            return targets
            
        else:
            return [instance(self.assignments)]
    
    def __str__(self):
        return '@template-target %s %s %s' % (
            self.name,
            self.template,
            self.assignments
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...

## VARIABLES LIST AND SUCH...
class List:
    def __init__(self, name, elements):
        self.name = name
        self.elements = elements

    def __str__(self):
        return '@list %s [%s]' % (
            self.name,
            ' '.join(self.elements)
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...

class Glob(List):
    def __init__(self, name, glob_pattern, elements):
        List.__init__(self, name, elements)
        self.glob_pattern = glob_pattern

    def __str__(self):
        return '@glob %s %s [%s]' % (
            self.name,
            self.glob_pattern,
            ' '.join(self.elements)
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...

class Shell(List):
    def __init__(self, name, shell_command, elements):
        List.__init__(self, name, elements)
        self.shell_command = shell_command

    def __str__(self):
        return '@shell %s %s [%s]' % (
            self.name,
            self.shell_command,
            ' '.join(self.elements)
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...

class Transform(List):
    def __init__(self, name, match_pattern, subs_pattern, input_list, elements):
        List.__init__(self, name, elements)
        self.match_pattern = match_pattern
        self.subs_pattern = subs_pattern
        self.input_list = input_list

    def __str__(self):
        return '@transform %s %s %s %s [%s]' % (
            self.name,
            self.match_pattern, self.subs_pattern,
            self.input_list,
            ' '.join(self.elements)
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...


## TASKS (TARGETS, FILES, AND ANYTHING THE WORKFLOW ACTUALLY SEES)
class Task:
    '''Abstract class for items in the workflow.'''
    def __init__(self, name, dependencies, wd):
        self.name = name
        self.dependencies = dependencies
        self.working_dir = wd
        self.is_dummy = False

    @property
    def should_run(self):
        '''Test if this task needs to be run. Used when scheduling a
        dependency graph, but the specifics of when and how a task should
        run depends on the subclass'''
        raise NotImplementedError()

    @property
    def can_execute(self):
        '''Flag used to indicate that a task can be executed.
        
        Sub-classes needs to implement it if they are capable of execution.
        By default the answer is false, so if a task that cannot execute
        is scheduled, we can report an error.'''
        return False
        
    @property
    def dummy(self):
        return self.is_dummy
        
    @property
    def execution_error(self):
        '''Should return an error message if a task that cannot execute
        is scheduled for execution.'''
        return 'Unknown task type cannot be executed.'
        
    @property
    def graphviz_shape(self):
        '''This is simply for UI. It provides a way of having different
        shapes for tasks when they are written to graphviz for workflow
        visualization.'''
        return 'box'


class SystemFile(Task):
    '''Class handling files that should be present on the system, i.e.
    files that are not generated by any targets but is specified as an
    input of one or more targets.'''
    
    def __init__(self, filename, wd):
        Task.__init__(self, filename, [], wd)
        self.is_dummy = True

    @property
    def file_exists(self):
        '''Check if the file exists. It is usually considered a major
        error if it doesn't since no target generates it.'''
        return _file_exists(_make_absolute_path(self.working_dir, self.name))

    @property
    def should_run(self):
        '''We should never actually run a system file, but we say yes when
        the file is missing so this is displayed in output.'''
        return not self.file_exists

    @property
    def execution_error(self):
        return 'The file "%s" is not generated by any target, ' \
            'so it must be present on the file server.' % self.name
            
    @property
    def graphviz_shape(self):
        return 'note'


class ExecutableTask(Task):
    '''Tasks that can be executed must provide this interface.'''

    def can_execute(self):
        return True
    
    @property
    def script_name(self):
        '''Where is the script for executing the task located?'''
        raise NotImplementedError()

    @property
    def job_name(self):
        '''Where is job-id file located if the task is executing?'''
        raise NotImplementedError()
    
    def write_script(self):
        '''Write the script for executing the task to disk.'''
        raise NotImplementedError()

    @property
    def graphviz_shape(self):
        return 'octagon'
        

class Target(ExecutableTask):
    '''Class handling targets. Stores the info for executing them.'''
    
    def __init__(self, name, input, output, pbs_options, flags, code, wd):
        # passing None as dependencies, 'cause Workflow will fill it in
        Task.__init__(self, name, None, wd)
        self.input = input
        self.output = output
        self.pbs_options = pbs_options
        self.flags = flags
        self.code = code
        
        if 'dummy' in self.flags and len(self.output) > 0:
            print 'Target %s is marked as a dummy target but has output files.'
            print 'Dummy targets will never be run so cannot produce output!'
            sys.exit(2)
        self.is_dummy = 'dummy' in self.flags
        

    @property
    def should_run(self):
        '''Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. '''
        
        
        if len(self.output) == 0:
            self.reason_to_run = \
                'Sinks (targets without output) should always run'
            return True # If we don't provide output, assume we always
                        # need to run.

               
        for outf in self.output:
            if not _file_exists(_make_absolute_path(self.working_dir, outf)):
                self.reason_to_run = \
                    'Output file %s is missing' % outf
                return True

        for inf in self.input:
            if not _file_exists(_make_absolute_path(self.working_dir,inf)):
                self.reason_to_run = \
                    'Input file %s is missing' % outf
                return True

        # If no file is missing, it comes down to the time stamps. If we
        # only have output and no input, we assume the output is up to
        # date. Touching files and adding input can fix this behaviour
        # from the user side but if we have a program that just creates
        # files we don't want to run it whenever someone needs that
        # output just because we don't have time stamped input.

        if len(self.input) == 0:
            self.reason_to_run = "We shouldn't run"
            return False

        # if we have both input and output files, check time stamps
        
        youngest_in_timestamp = None
        youngest_in_filename = None
        for inf in self.input:
            timestamp = _get_file_timestamp(_make_absolute_path(self.working_dir,inf))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_filename = inf
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        oldest_out_timestamp = None
        oldest_out_filename = None
        for outf in self.output:
            timestamp = _get_file_timestamp(_make_absolute_path(self.working_dir,outf))
            if oldest_out_timestamp is None \
                    or oldest_out_timestamp > timestamp:
                oldest_out_filename = outf
                oldest_out_timestamp = timestamp
        assert oldest_out_timestamp is not None
        
        # The youngest in should be older than the oldest out
        if youngest_in_timestamp >= oldest_out_timestamp:
            # we have a younger in file than an outfile
            self.reason_to_run = 'Infile %s is younger than outfile %s' %\
                (youngest_in_filename, oldest_out_filename)
            return True
        else:
            self.reason_to_run = 'Youngest infile %s is older than '\
                                 'the oldest outfile %s' % \
                (youngest_in_filename, oldest_out_filename)
            return False    
            
        assert False, "We shouldn't get here"
        

    @property
    def script_dir(self):
        return _make_absolute_path(self.working_dir,'.scripts')
    
    @property
    def jobs_dir(self):
        return _make_absolute_path(self.working_dir,'.jobs')

    @property
    def script_name(self):
        # Escape name to make a file name...
        escaped_name = _escape_file_name(self.name)
        return _make_absolute_path(self.script_dir, escaped_name)

    @property
    def job_name(self):
        # Escape name to make a file name...
        escaped_name = _escape_file_name(self.name)
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
        for options in self.pbs_options:
            print >> f, '#PBS', options
        print >> f
        
        print >> f, '# GWF generated code ...'
        print >> f, 'cd %s' % self.working_dir
        print >> f

        print >> f, '# Script from workflow'
        print >> f, self.code

    @property
    def job_in_queue(self):
        if not _file_exists(self.job_name):
            return False
        else:
            # FIXME: check status
            jobid = self.jobID
            stat = subprocess.Popen(['qstat','-f',jobid],
                                    stdout=subprocess.PIPE)
            for line in stat.stdout:
                line = line.strip()
                print line
                if line.startswith('job_state'):
                    return True
            return False
    
    @property
    def jobID(self):
        if _file_exists(self.job_name):
            return open(self.job_name).read().strip()
        else:
            return None
        

    @property
    def graphviz_shape(self):
        if len(self.input) == 0 and len(self.output) > 0:
            return 'invhouse'
        if len(self.output) == 0 and len(self.input) > 0:
            return 'house'
        
        return 'octagon'

    def __str__(self):
        return '@target %s, input(%s) -> output(%s)' % (
            self.name,
            ' '.join(self.input),
            ' '.join(self.output)
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy 
    				   # for printing output when testing...


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
                    print "Transformation list %s uses input list %s the doesn't exist."%\
                        (cmd.name, input_list_name)
                    sys.exit(2)
                
                input_list = self.lists[input_list_name]
                
                cmd.elements = [re.sub(cmd.match_pattern,cmd.subs_pattern,input)
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
                    print 'Instantiated template %s has the same name as an existing target' %\
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
            target.input  = expand_lists(target.input)
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
                import sys ; sys.exit(2)
        
            # make sure we have the scripts for the jobs we want to
            # execute!
            job.task.write_script()

            dependent_tasks = set(node.name
                                  for _,node in job.dependencies
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
                'echo', ('$%s'%job.name), '>', job.task.job_name]))

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
                import sys ; sys.exit(2)
        
            # make sure we have the scripts for the jobs we want to
            # execute!
            job.task.write_script()

            script = open(job.task.script_name, 'r').read()
            script_commands.append('# computing %s' % job.name)
            script_commands.append(script)

        return '\n'.join(script_commands)


