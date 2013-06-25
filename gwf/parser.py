'''Functionality for parsing a workflow specification file.'''

import sys
import os.path
import re
from workflow import Template, TemplateTarget, Target, Workflow

# working dir is provided to the parse in case we want to allow
# recursive workflows or something in the future... and I need it for
# the files in Targets() at some point, and I might as well provide it
# here rather than setting it in the workflow class later...

def parse_target(target_code, working_dir):
    lines = target_code.split('\n')
    _, name = lines[0].split()

    input = []
    output = []
    pbs = []

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

            else:
                print "Unknown opcode %s in target %s." % (line, name)
                sys.exit(2)

        else:
            # assumes everything after the opcodes is the script
            break
                
    code = '\n'.join(lines[i:])
            
    return Target(name, input, output, pbs, code, working_dir)

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
    header, fluff = code.split('\n',1)
    header_objects = header.split()
    name = header_objects[1]
    template = header_objects[2]
    parameter_assignments = header_objects[3:]
    assignments = dict()
    for assignment in parameter_assignments:
        key, val = assignment.split('=')
        assignments[key.strip()] = val.strip()
        
    return TemplateTarget(name, working_dir, template, assignments)


PARSERS = {'target': parse_target,
           'template': parse_template,
           'template-target': parse_template_target}

def parse(fname):
    '''Parse up the workflow in "fname".'''

    working_dir = os.path.dirname(os.path.realpath(fname))
    workflow_text = open(fname).read()
    commands = workflow_text.split('\n@')[1:]

    parsed_commands = []
    for cmd in commands:
        opcode,_ = re.split('\s',cmd,1)
        
        if opcode not in PARSERS:
            print 'Unknown task type @%s.' % opcode
            sys.exit(2)

        parsed_commands.append(PARSERS[opcode](cmd, working_dir))

    templates = dict()
    targets = dict()
    template_targets = dict()
    for cmd in parsed_commands:
        # FIXME: probably shouldn't hardwire a test for the type of task here!
        
        if isinstance(cmd, Template):
            if cmd.name in templates:
                print 'Template %s appears more than once in the workflow.' % cmd.name
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
                print 'Task %s appears more than once in the workflow.' % cmd.name
                sys.exit(2)
            
            targets[cmd.name] = cmd

    return Workflow(templates, targets, template_targets, working_dir)

