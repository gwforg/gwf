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

def _format_options(options):
    # FIXME: Check if the options are valid!
    options_list = []
    for key, value in options.items():
        # FIXME: format the value based on type
        if hasattr(value, '__iter__'):
            # Any list or such we assume should just be space separated strings
            value = ' '.join(map(str, value))
        
        if type(value) == bool:
            # booleans should be set as flags if True and ignored otherwise
            if value:
                value = ''  # becomes the empty string in the ":foo" syntax
            else:
                continue    # no :foo set at all
            
        options_list.append(':{key} {value}'.format(key=key, value=value))
    return '\n'.join(options_list)

class target(object):
    def __init__(self, name, **options):
        self.name = name
        self.options = options
        
    def __lshift__(self, spec):
        ## FIXME: Should build a node for the workflow graph instead of outputting a textual target
        
        options = self.options
        if isinstance(spec, _TemplateInstance):
            spec = spec.spec
            options = dict(options.items() + spec.options.items())
    
        if self.name in gwf_workflow.ALL_TARGETS:
            print 'Warning: Target', self.name, 'defined more than once.'
        new_target = gwf_workflow.Target(self.name, options, spec)
        gwf_workflow.ALL_TARGETS[self.name] = new_target
    
        #namespace['_formatted_options_'] = _format_options(options)
        #print '@target {_target_name_}\n{_formatted_options_}\n{_spec_}\n'.format(**namespace)

