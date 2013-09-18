import os
from mako.template import Template

GWF_IMPORTS = '''
import os
import re
import glob as _glob
import subprocess
'''

SHELL_FUNCTIONS = '''
def glob(pattern):
    if pattern.startswith('/'):
        glob_pattern = pattern
    else:
        glob_pattern = os.path.join(working_dir,pattern)
    return _glob.glob(glob_pattern)

def shell(command):
    return subprocess.check_output(command, shell=True).split()

'''

def preprocess(fname):
    working_dir = os.path.dirname(os.path.realpath(fname))
    GWF_FORMATTING_IMPORTS = [
        GWF_IMPORTS, 
        'working_dir = "%s"' % working_dir,
        SHELL_FUNCTIONS,
        ]
    template = Template(filename=fname, imports=GWF_FORMATTING_IMPORTS)
    workflow_text = template.render()
    return workflow_text