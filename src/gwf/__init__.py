import os
import os.path
import sys
import inspect
import glob as _glob
import subprocess

import gwf_workflow

## Useful helper functions...


def glob(pattern):
    if pattern.startswith('/'):
        glob_pattern = pattern
    else:
        filename = inspect.getfile(sys._getframe(1))
        working_dir = os.path.dirname(os.path.realpath(filename))
        glob_pattern = os.path.join(working_dir, pattern)
    return _glob.glob(glob_pattern)


def shell(command):
    return subprocess.check_output(command, shell=True).split()


## Templates
class template(object):
    def __init__(self, **options):
        self.spec = None
        self.options = options

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __call__(self, **substitutions):
        def substitute(s):
            if type(s) == str:
                return s.format(**substitutions)
            elif hasattr(s, '__iter__'):
                return [substitute(x) for x in s]
            else:
                return s

        formatted_options = [(key, substitute(val)) for key, val in self.options.items()]
        options = dict(formatted_options)
        return options, self.spec.format(**substitutions)


# Targets...
class target(object):
    def __init__(self, name, **options):
        self.name = name
        self.options = options

    def __lshift__(self, spec):
        options = self.options

        if isinstance(spec, tuple):
            options = dict(options.items() + spec[0].items())
            spec = spec[1]

        if self.name in gwf_workflow.ALL_TARGETS:
            print 'Warning: Target', self.name, 'defined more than once.'

        new_target = gwf_workflow.Target(self.name, options, spec)
        gwf_workflow.ALL_TARGETS[self.name] = new_target


# This is a bit of a hack for now, but I need a way of specifying the backend
GWF_BACKEND="PSB"
def set_backend(backend):
    global GWF_BACKEND
    GWF_BACKEND = backend
