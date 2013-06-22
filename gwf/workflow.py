'''Classes representing a workflow.'''

# FIXME: handle errors

class Target:
    '''Class handling targets. Stores the info for executing them.'''
    
    def __init__(self, name, input, output, code):
        self.name = name
        self.input = input
        self.output = output
        self.code = code

        # dependencies are the targets that provide files.
        self.dependencies = None # will be filled in by Workflow

        # system_files are expected to be existing, since
        # we don't have targets providing them.
        self.system_files = None # will be filled in by Workflow

    def __str__(self):
        sf = ''
        dep = ''
        if self.system_files:
            sf = '%s must be present on the file server' % \
                ', '.join(self.system_files)
        if self.dependencies:
            dep = ','.join([('target %s should provide %s ' % (tg.name,of))
                            for of,tg in self.dependencies])
        return '@target %s, input(%s) -> output(%s), %s [%s; %s]' % (
            self.name,
            ' '.join(self.input),
            ' '.join(self.output),
            '%d lines of code' % len(self.code.split('\n')),
            sf, dep
            )
    __repr__ = __str__


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
