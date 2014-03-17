import os
import os.path
import re
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
        fname = inspect.getfile(sys._getframe(1))
        working_dir = os.path.dirname(os.path.realpath(fname))
        glob_pattern = os.path.join(working_dir,pattern)
    return _glob.glob(glob_pattern)

def shell(command):
    return subprocess.check_output(command, shell=True).split()


## Templates
class _TemplateInstance(object):
    def __init__(self, spec, options):
        self.spec = spec
        self.options = options

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
        
        formatted_options = [(key,substitute(val)) for key,val in self.options.items()]
        options = dict(self.options.items() + formatted_options)
        return _TemplateInstance(self.spec.format(**substitutions), options)

     
# Targets...
class target(object):
    def __init__(self, name, **options):
        self.name = name
        self.options = options
        
    def __lshift__(self, spec):
        options = self.options
        if isinstance(spec, _TemplateInstance):
            spec = spec.spec
            options = dict(options.items() + spec.options.items())
    
        if self.name in gwf_workflow.ALL_TARGETS:
            print 'Warning: Target', self.name, 'defined more than once.'
        new_target = gwf_workflow.Target(self.name, options, spec)
        gwf_workflow.ALL_TARGETS[self.name] = new_target
    
