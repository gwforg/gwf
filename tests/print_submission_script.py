from gwf.parser import parse

if __name__ == '__main__':
    # FIXME: make a better test...
    import sys
    workflow = parse(sys.argv[1])
    script = workflow.get_submission_script(sys.argv[2])

    print '# To submit target %s the following script will be called.' %\
        sys.argv[2]
    print script


