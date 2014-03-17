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

        self.input = []
        self.output = []
        self.pbs = []

        known_options = ('input', 'output', 'pbs')
        for k in options.keys():
            if k not in known_options:
                print 'Warning:, Target', self.name, 'has unknown option', k

        if 'input' in options:
            self.input = _list(options['input'])
        if 'output' in options:
            self.output = _list(options['output'])
        if 'pbs' in options:
            self.pbs = _list(options['pbs'])

    def __str__(self):
        return '''@target {name}\n:input {input}\n:output {output}\n:pbs {pbs}\n\n{spec}'''.format(
            name=self.name,
            input=' '.join(self.input), output=' '.join(self.output),
            pbs=' '.join(self.pbs),
            spec=self.spec
        )

    __repr__ = __str__


ALL_TARGETS = {}

