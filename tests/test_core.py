import os.path
import unittest
from unittest.mock import patch

from gwf import Workflow, Target
from gwf.exceptions import GWFException


class TestWorkflow(unittest.TestCase):

    def setUp(self):
        self.workflow = Workflow()

    def test_adding_a_target_makes_it_available_to_the_workflow(self):
        self.workflow.target('TestTarget', inputs=[], outputs=[], spec='')
        self.assertIn('TestTarget', self.workflow.targets)

    def test_adding_two_targets_with_the_same_names_should_raise_an_exception(self):
        self.workflow.target('TestTarget', inputs=[], outputs=[], spec='')
        with self.assertRaises(GWFException):
            self.workflow.target('TestTarget', inputs=[], outputs=[], spec='')


class TestTarget(unittest.TestCase):

    def create_test_target(name='TestTarget',
                           inputs=[],
                           outputs=[],
                           options={},
                           working_dir=''):
        return Target(name, inputs, outputs, options, working_dir)

    def test_relative_input_paths_are_normalized(self):
        target = self.create_test_target(
            inputs=['test_input1.txt', 'test_input2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(os.path.isabs(target.inputs[0]))
        self.assertTrue(os.path.isabs(target.inputs[1]))

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/some/path'))

    def test_relative_output_paths_are_normalized(self):
        target = self.create_test_target(
            outputs=['test_output1.txt', 'test_output2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(os.path.isabs(target.outputs[0]))
        self.assertTrue(os.path.isabs(target.outputs[1]))

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/some/path'))

    def test_absolute_input_paths_are_not_normalized(self):
        target = self.create_test_target(
            inputs=['test_input1.txt', '/other/path/test_input2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_absolute_output_paths_are_not_normalized(self):
        target = self.create_test_target(
            inputs=['test_output1.txt', '/other/path/test_output2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_target_without_outputs_is_a_sink(self):
        target = self.create_test_target()
        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = self.create_test_target(
            outputs=['test_output1.txt', 'test_output2.txt']
        )
        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = self.create_test_target()
        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = self.create_test_target(
            inputs=['test_input1.txt', 'test_input2.txt']
        )
        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = self.create_test_target() << 'this is a spec'
        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, 'this is a spec')


class TestDependencies(unittest.TestCase):
    pass


class TestGetExecutionSchedule(unittest.TestCase):
    pass
