'''Functionality for parsing a workflow specification file.'''

import re
from workflow import Target, Workflow

def parse_target(target_code):
    lines = target_code.split('\n')
    _, name = lines[0].split()

    input = []
    output = []

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

            else:
                assert False, 'Unknown opcode %s' % line
                # FIXME: Proper error handling
        else:
            # assumes everything after the opcodes is the script
            code = '\n'.join(lines[i:])
            break
                
    return Target(name, input, output, code)

PARSERS = {'target': parse_target}

def parse(fname):
    '''Parse up the workflow in "fname".'''
    # FIXME: add error handling!
    workflow_text = open(fname).read()
    commands = workflow_text.split('\n@')[1:]

    workflow = []
    for cmd in commands:
        opcode,_ = re.split('\s',cmd,1)
        assert opcode in PARSERS # FIXME: add error handling

        workflow.append(PARSERS[opcode](cmd))

    targets = dict()
    for cmd in workflow:
        if isinstance(cmd, Target):
            # FIXME: error handling
            assert cmd.name not in targets
            targets[cmd.name] = cmd

    return Workflow(targets)

if __name__ == '__main__':
    # FIXME: make a better test...
    import sys
    workflow = parse(sys.argv[1])

    print 'Workflow contains %d targets:' % len(workflow.targets)
    for name,target in workflow.targets.items():
        print target
    print
