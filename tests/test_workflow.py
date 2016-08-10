import unittest

from gwf import Workflow, Target


class TestWorkflow(unittest.TestCase):
    pass


class TestTarget(unittest.TestCase):

    def test_input_paths_are_normalized(self):
        self.assertTrue(False)

    def test_output_paths_are_normalized(self):
        self.assertTrue(False)

    def test_target_without_outputs_is_a_sink(self):
        self.assertTrue(False)

    def test_target_without_inputs_is_a_source(self):
        self.assertTrue(False)


class TestDependencies(unittest.TestCase):
    pass


class TestGetExecutionSchedule(unittest.TestCase):
    pass
