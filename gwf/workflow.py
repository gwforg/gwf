'''Classes representing a workflow.'''

# FIXME: handle errors

import os.path, time

def _file_exists(fname):
    return os.path.exists(fname)

def _get_file_timestamp(fname):
    return time.ctime(os.path.getmtime(fname))

class Target:
    '''Class handling targets. Stores the info for executing them.'''
    
    def __init__(self, name, input, output, code, wd):
        self.name = name
        self.input = input
        self.output = output
        self.code = code
        self.working_dir = wd

        # dependencies are the targets that provide files.
        self.dependencies = None # will be filled in by Workflow

        # system_files are expected to be existing, since
        # we don't have targets providing them.
        self.system_files = None # will be filled in by Workflow

    def check_system_files(self):
        '''Make sure that the files we don't know how to make are
        actually on the file server.'''

        for sysfile in self.system_files:
           assert _file_exists(self.working_dir+'/'+sysfile), \
           "There is no target providing %s. It must be present!" % sysfile

    def should_run(self):
        '''Test if this target needs to be run.'''
        
        if len(self.output) == 0:
            return True # If we don't provide output, assume we always
                        # need to run.
        if len(self.input) == 0:
            return True # If we don't have any input, we also assume
                        # we have to run

        # If we have both input and output, we need to run if any of
        # the input files are newer than any of the output
        # files... and if any of the output or inputs are missing, we
        # definitely have to run.
        
        for outf in self.output:
            # FIXME: Handle absolute paths...
            if not _file_exists(self.working_dir+'/'+outf):
                return True

        for inf in self.input:
            if not _file_exists(self.working_dir+'/'+inf):
                # FIXME: error handling
                assert inf not in self.system_files, \
                    "There is no target providing %s. It must be present!" % inf
                return True

        # if we have all input and output files, check time stamps
        in_timestamp = min(_get_file_timestamp(self.working_dir+'/'+inf)
                           for inf in self.input)
        out_timestamp = min(_get_file_timestamp(self.working_dir+'/'+outf)
                            for outf in self.output)

        return in_timestamp > out_timestamp


    def __str__(self):
        sf = ''
        dep = ''
        if self.system_files:
            sf = '%s must be present on the file server' % \
                ', '.join(self.system_files)
        if self.dependencies:
            dep = ','.join([('target %s should provide %s ' % (tg.name,of))
                            for of,tg in self.dependencies])
        return '@target %s, input(%s) -> output(%s), %s [%s; %s] %s' % (
            self.name,
            ' '.join(self.input),
            ' '.join(self.output),
            '%d lines of code' % len(self.code.split('\n')),
            sf, dep,
            'Should run' if self.should_run() else "Doesn't need to run"
            )
    __repr__ = __str__ # not really the correct use of __repr__ but easy for printing output when testing...


class Workflow:
    '''Class representing a workflow.'''

    def __init__(self, targets):
        self.targets = targets

        # collect the output files so we know who can build them.
        self.providers = dict()
        for target in self.targets.values():
            for output_file in target.output:
                assert output_file not in self.providers
                self.providers[output_file] = target

        # now get dependencies for each target...
        for target in self.targets.values():
            dependencies = []
            system_files = []
            for input_file in target.input:
                if input_file in self.providers:
                    dependencies.append((input_file,
                                         self.providers[input_file]))
                else:
                    system_files.append(input_file)
            target.dependencies = dependencies
            target.system_files = system_files
