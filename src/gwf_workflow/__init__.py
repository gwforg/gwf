'''
Module for building and executing workflows.
'''

class Target(object):
    def __init__(self, name, options, spec):
        self.name = name
        self.options = options
        self.spec = spec

ALL_TARGETS = {}

