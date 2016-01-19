"""
Module for building and executing work flows.
"""

import sys
import inspect
import os.path


# Helper function for wrapping singletons in lists
def _list(x):
    """Wrap x as a singleton in a list if it isn't a list already."""
    if hasattr(x, '__iter__'):
        return list(x)
    else:
        return [x]


# Internal representation of targets.
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

        known_options_without_defaults = set(['queue', 'account'])

        for k in options.keys():
            if k in self.options or k in known_options_without_defaults:
                self.options[k] = options[k]
            else:
                print 'Warning:, Target', self.name, 'has unknown option', k

        # handle that input and output can be both lists and single file names
        self.options['input'] = filter(None, _list(self.options['input'])) # safer in case someone writes input='' meaning "no input"
        self.options['output'] = filter(None, _list(self.options['output']))


# This global variable will hold all the targets after a workflow script has completed.
# gwf will use this list for its further processing.
ALL_TARGETS = {}


# This will be set in the gwf script and refer to the grid backend used.
BACKEND = None