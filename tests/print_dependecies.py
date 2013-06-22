from gwf.parser import parse

if __name__ == '__main__':
    # FIXME: make a better test...
    import sys
    workflow = parse(sys.argv[1])
    target = workflow.targets[sys.argv[2]]

    target.get_dependencies().print_dependency_graph()
