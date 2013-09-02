'''Functionality for parsing a workflow specification file.'''

import sys
import os.path
import re
import glob
import subprocess
from workflow import List
from workflow import Glob, Shell, Transform, SubGraph
from workflow import Template, TemplateTarget
from workflow import Target
from workflow import Workflow

# working dir is provided to the parse in case we want to allow
# recursive workflows or something in the future... and I need it for
# the files in Targets() at some point, and I might as well provide it
# here rather than setting it in the workflow class later...

def parse_comment(code, working_dir):
    pass

def parse_target(target_code, working_dir):
    lines = target_code.split('\n')
    _, name = lines[0].split()

    input = []
    output = []
    pbs = []
    flags = []

    for i in xrange(1,len(lines)):
        line = lines[i]
        if line.startswith(':'):
            # the line is an opcode for targets
            if line.startswith(':input'):
                files = line.split()[1:]
                input.extend(files)

            elif line.startswith(':output'):
                files = line.split()[1:]
                output.extend(files)
                
            elif line.startswith(':pbs'):
                options = line.split()[1:]
                pbs.append(' '.join(options))

            elif line.startswith(':dummy'):
                flags.append('dummy')

            else:
                print "Unknown opcode %s in target %s." % (line, name)
                sys.exit(2)

        else:
            # assumes everything after the opcodes is the script
            break
                
    code = '\n'.join(lines[i:])
            
    return Target(name, input, output, pbs, flags, code, working_dir)

def parse_template(template_code, working_dir):
    header, code = template_code.split('\n',1)
    header_objects = header.split()
    name = header_objects[1]
    parameters = header_objects[2:]
    if len(parameters) == 0:
        print 'Template %s has no parameters. That must be an error.' % name
        sys.exit(2)
        
    return Template(name, working_dir, parameters, code)

def parse_template_target(code, working_dir):
    elements = code.split()
    name = elements[1]
    template = elements[2]
    parameter_assignments = elements[3:]
    assignments = dict()
    for assignment in parameter_assignments:
        key, val = assignment.split('=')
        assignments[key.strip()] = val.strip()
        
    return TemplateTarget(name, working_dir, template, assignments)


def parse_list(code, working_dir):
    elements = code.split()
    name = elements[1]
    values = elements[2:]
    return List(name, values)
    
def parse_glob(code, working_dir):
    elements = code.split()
    if len(elements) != 3:
        print 'Malformed glob "%s"' % code
        sys.exit(2)
    
    name = elements[1]
    pattern = elements[2]
    if pattern.startswith('/'):
        glob_pattern = pattern
    else:
        glob_pattern = os.path.join(working_dir,pattern)
    elements = glob.glob(glob_pattern)
    
    return Glob(name, glob_pattern, elements)

def parse_shell(code, working_dir):
    _,name,command = re.split('\s+',code,2)
    command = command.strip()
    shell_result = subprocess.check_output(command, shell=True).split()
    return Shell(name, command, shell_result)
    
def parse_transform(code, working_dir):
    elements = code.split()
    if len(elements) != 5:
        print 'Malformed transform "%s"' % code
        sys.exit(2)
    
    name = elements[1]
    match_pattern = elements[2]
    subs_pattern = elements[3]
    input_list = elements[4]
    elements = None # we can't transform yet, but the workflow will do it soon!
    
    return Transform(name, match_pattern, subs_pattern, input_list, elements)

def parse_subgraph(code, working_dir):
    elements = code.split()
    if len(elements) != 2:
        print 'Malformed subgraph "%s"' % code
        sys.exit(2)
    return SubGraph(elements[1].strip())
    
    
PARSERS = {'target': parse_target,
           'template': parse_template,
           'template-target': parse_template_target,
           'list': parse_list, 
           'glob': parse_glob,
           'shell': parse_shell,
           'transform': parse_transform,
           'comment': parse_comment,
           'subgraph': parse_subgraph,
            }

def command_iter(fname):
    '''Iterator for the commands in a file.'''
    working_dir = os.path.dirname(os.path.realpath(fname))
    try:
        workflow_text = open(fname).read()
    except:
        print "Problems opening file '%s'." % fname
        sys.exit(2)
        
    commands = workflow_text.split('\n@')[1:]
    for cmd in commands:
        opcode,_ = re.split('\s',cmd,1)
        
        if opcode not in PARSERS:
            print 'Unknown task type @%s.' % opcode
            sys.exit(2)

        next_parsed_cmd = PARSERS[opcode](cmd, working_dir)
        if isinstance(next_parsed_cmd, SubGraph):
            for sub_cmd in command_iter(os.path.join(working_dir, next_parsed_cmd.filename)):
                yield sub_cmd
        else:
            yield next_parsed_cmd

def parse(fname):
    '''Parse up the workflow in "fname".'''
    working_dir = os.path.dirname(os.path.realpath(fname))

    lists = dict()
    templates = dict()
    targets = dict()
    template_targets = dict()
    for cmd in command_iter(fname):
        
        if isinstance(cmd, List):
            if cmd.name in lists:
                print 'List %s appears more than once in the workflow.' % \
                    cmd.name
                sys.exit(2)
            lists[cmd.name] = cmd
        
        if isinstance(cmd, Template):
            if cmd.name in templates:
                print 'Template %s appears more than once in the workflow.' % \
                    cmd.name
                sys.exit(2)
            templates[cmd.name] = cmd
        
        if isinstance(cmd, TemplateTarget):
            if cmd.name in template_targets:
                print 'Template target %s appears more than once in the worlflow.' % \
                    cmd.name
                sys.exit(2)
            template_targets[cmd.name] = cmd
        
        if isinstance(cmd, Target):

            if cmd.name in targets:
                print 'Task %s appears more than once in the workflow.' % \
                    cmd.name
                sys.exit(2)
            targets[cmd.name] = cmd

    return Workflow(lists, templates, targets, template_targets, working_dir)

