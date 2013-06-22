'''Classes representing a workflow.'''

class Target:
    '''Class handling targets. Stores the info for executing them.'''
    
    def __init__(self, name, input, output, code):
        self.name = name
        self.input = input
        self.output = output
        self.code = code


    def __str__(self):
        return '@target %s, input(%s) -> output(%s), %s' % (
            self.name,
            ' '.join(self.input),
            ' '.join(self.output),
            '%d lines of code' % len(self.code.split('\n'))
            )
    __repr__ = __str__


class Workflow:
    '''Class representing a workflow.'''

    def __init__(self, targets):
        self.targets = targets


        
