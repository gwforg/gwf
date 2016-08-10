import os.path
import unittest

from gwf import Workflow, Target


class TestWorkflow(unittest.TestCase):
    pass


class TestTarget(unittest.TestCase):

    def test_relative_input_paths_are_normalized(self):
        target = Target(name='TestTarget',
                        inputs=['test_input1.txt', 'test_input2.txt'],
                        outputs=[],
                        options={},
                        working_dir='/some/path')

        self.assertTrue(os.path.isabs(target.inputs[0]))
        self.assertTrue(os.path.isabs(target.inputs[1]))

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/some/path'))

    def test_relative_output_paths_are_normalized(self):
        target = Target(name='TestTarget',
                        inputs=[],
                        outputs=['test_output1.txt', 'test_output2.txt'],
                        options={},
                        working_dir='/some/path')

        self.assertTrue(os.path.isabs(target.outputs[0]))
        self.assertTrue(os.path.isabs(target.outputs[1]))

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/some/path'))

    def test_absolute_input_paths_are_not_normalized(self):
        self.assertTrue(False)

    def test_absolute_output_paths_are_not_normalized(self):
        self.assertTrue(False)

    def test_target_without_outputs_is_a_sink(self):
        target = Target(name='TestTarget',
                        inputs=[],
                        outputs=[],
                        options={},
                        working_dir='')

        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = Target(name='TestTarget',
                        inputs=[],
                        outputs=['test_output1.txt', 'test_output2.txt'],
                        options={},
                        working_dir='')

        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = Target(name='TestTarget',
                        inputs=[],
                        outputs=[],
                        options={},
                        working_dir='')

        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = Target(name='TestTarget',
                        inputs=['test_input1.txt', 'test_input2.txt'],
                        outputs=[],
                        options={},
                        working_dir='')

        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = Target(name='TestTarget',
                        inputs=[],
                        outputs=[],
                        options={},
                        working_dir='') << 'this is a spec'

        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, 'this is a spec')


class TestDependencies(unittest.TestCase):
    pass


class TestGetExecutionSchedule(unittest.TestCase):
    pass
