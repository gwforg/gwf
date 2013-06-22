from gwf.parser import parse

if __name__ == '__main__':
    # FIXME: make a better test...
    import sys
    workflow = parse(sys.argv[1])
    target = workflow.targets[sys.argv[2]]

    schedule = target.get_dependencies().schedule()
    print schedule
    for node in schedule:
        node.target.execute_locally()

