import unittest
from unittest.mock import patch

from gwf import Target, Workflow, inputs, outputs


class GWFTestCase(unittest.TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing


def create_test_target(name='TestTarget',
                       input_files=None,
                       output_files=None,
                       options=None,
                       workflow=None):
    """A factory for `Target` objects."""
    if input_files is None:
        input_files = []
    if output_files is None:
        output_files = []
    if options is None:
        options = {}
    if workflow is None:
        workflow = Workflow(working_dir='/some/path')
    target = Target(name, options, workflow)
    target = target << inputs(*input_files) << outputs(*output_files)
    return target
