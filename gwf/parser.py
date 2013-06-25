'''Functionality for parsing a workflow specification file.'''

import sys
import os.path
import re
from workflow import Target, Workflow

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

PARSERS = {'target': parse_target}

def parse(fname):
    '''Parse up the workflow in "fname".'''

    working_dir = os.path.dirname(os.path.realpath(fname))
    workflow_text = open(fname).read()
    commands = workflow_text.split('\n@')[1:]

    workflow = []
    for cmd in commands:
        opcode,_ = re.split('\s',cmd,1)
        
        if opcode not in PARSERS:
            print 'Unknown task type @%s.' % opcode
            sys.exit(2)

        workflow.append(PARSERS[opcode](cmd, working_dir))

    targets = dict()
    for cmd in workflow:
        # FIXME: probably shouldn't hardwire a test for the type of task here!
        if isinstance(cmd, Target):

            if cmd.name in targets:
                print 'Task %s appears more than once in the workflow.' % cmd.name
                sys.exit(2)
            
            targets[cmd.name] = cmd

    return Workflow(targets, working_dir)

