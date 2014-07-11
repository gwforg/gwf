"""
Module for building and executing work flows.
"""

import sys
import inspect
import os.path


def _list(x):
    """Wrap x as a singleton in a list if it isn't a list already."""
    if hasattr(x, '__iter__'):
        return list(x)
    else:
        return [x]


class Target(object):
    def __init__(self, name, options, spec):
        self.name = name
        self.spec = spec

        filename = inspect.getfile(sys._getframe(2))
        self.working_dir = os.path.dirname(os.path.realpath(filename))

        self.options = {
            'input': [],
            'output': [],
            'nodes': 1,
            'cores': 1,
            'memory': "4g",
            'walltime': '120:00:00',
        }

        for k in options.keys():
            if k in self.options:
                self.options[k] = options[k]
            else:
                print 'Warning:, Target', self.name, 'has unknown option', k

        # handle that input and output can be both lists and single file names
        self.options['input'] = _list(self.options['input'])
        self.options['output'] = _list(self.options['output'])

    def __str__(self):
        return '''@target {name}\n:input {input}\n:output {output}\n:pbs {pbs}\n\n{spec}'''.format(
            name=self.name,
            input=' '.join(self.input), output=' '.join(self.output),
            pbs=' '.join(self.pbs),
            spec=self.spec
        )

    __repr__ = __str__


ALL_TARGETS = {}

