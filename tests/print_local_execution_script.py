from gwf.parser import parse

if __name__ == '__main__':
    # FIXME: make a better test...
    import sys
    workflow = parse(sys.argv[1])
    target = workflow.targets[sys.argv[2]]

    schedule = target.get_dependencies().schedule()
    print '# To execute target %s the following script will be called.' %\
        target.name
    for n in schedule:
        print n.target.local_execution_script()
    print

