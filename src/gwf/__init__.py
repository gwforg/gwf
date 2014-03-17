import os
import os.path
import re
import sys
import inspect
import glob as _glob
import subprocess

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
        
    def __lshift__(self, _spec_):
        ## FIXME: Should build a node for the workflow graph instead of outputting a textual target
        
        namespace = dict()
        namespace['_target_name_'] = self.name
        namespace['_spec_'] = _spec_
        options = self.options
        
        if isinstance(_spec_, _TemplateInstance):
            namespace['_spec_'] = _spec_.spec
            options = dict(options.items() + _spec_.options.items())
    
        namespace['_formatted_options_'] = _format_options(options)
                                                                         
        print '@target {_target_name_}\n{_formatted_options_}\n{_spec_}\n'.format(**namespace)

