import unittest
from unittest.mock import patch

from gwf import Target, Workflow


class GWFTestCase(unittest.TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing


def create_test_target(name='TestTarget', inputs=None, outputs=None, options=None, workflow=None):
    """A factory for `Target` objects."""
    if inputs is None:
        inputs = []
    if outputs is None:
        outputs = []
    if options is None:
        options = {}
    if workflow is None:
        workflow = Workflow(working_dir='/some/path')
    return Target(name, inputs, outputs, options, workflow)
