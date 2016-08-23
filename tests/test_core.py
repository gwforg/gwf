import os.path
import unittest
from unittest.mock import Mock, create_autospec, patch

from gwf.core import PreparedWorkflow, Target, Workflow
from gwf.exceptions import (CircularDependencyError,
                            FileProvidedByMultipleTargetsError,
                            FileRequiredButNotProvidedError, TargetExistsError)


def create_test_target(name='TestTarget', inputs=[], outputs=[], options={}, working_dir=''):
    """A factory for `Target` objects."""
    return Target(name, inputs, outputs, options, working_dir)


class TestWorkflow(unittest.TestCase):

    def test_target_with_no_input_has_empty_inputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget')
        self.assertListEqual(target.inputs, [])

    def test_target_with_no_output_has_empty_outputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget')
        self.assertListEqual(target.outputs, [])

    def test_adding_a_target_makes_it_available_to_the_workflow(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=[], spec='')

        self.assertIn('TestTarget', workflow.targets)

    def test_adding_two_targets_with_the_same_names_should_raise_an_exception(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=[], spec='')

        with self.assertRaises(TargetExistsError):
            workflow.target('TestTarget', inputs=[], outputs=[], spec='')

    def test_including_workflow_object_should_extend_including_workflow(self):
        workflow = Workflow()
        workflow.target('TestTarget1', inputs=[], outputs=[])

        other_workflow = Workflow()
        other_workflow.target('TestTarget2', inputs=[], outputs=[])
        other_workflow.target('TestTarget3', inputs=[], outputs=[])

        workflow.include_workflow(other_workflow)

        self.assertIn('TestTarget1', workflow.targets)
        self.assertIn('TestTarget2', workflow.targets)
        self.assertIn('TestTarget3', workflow.targets)

    def test_including_workflow_object_with_target_with_existing_name_should_raise_an_exception(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=[])

        other_workflow = Workflow()
        other_workflow.target('TestTarget', inputs=[], outputs=[])

        with self.assertRaises(TargetExistsError):
            workflow.include_workflow(other_workflow)

    @patch('gwf.core.import_object')
    def test_including_workflow_path_import_object_and_include_workflow_into_current_workflow(self, mock_import_object):
        workflow = Workflow()
        workflow.target('TestTarget1', inputs=[], outputs=[])

        other_workflow = Workflow()
        other_workflow.target('TestTarget2', inputs=[], outputs=[])
        other_workflow.target('TestTarget3', inputs=[], outputs=[])

        mock_import_object.return_value = other_workflow

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include_path('/path/to/other_workflow.py')

            mock_import_object.assert_called_once_with(
                '/path/to/other_workflow.py'
            )

            mock_include_workflow.assert_called_once_with(
                other_workflow
            )

    def test_including_workflow_instance_dispatches_to_include_workflow(self):
        workflow = Workflow()
        other_workflow = Workflow()

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include(other_workflow)
            mock_include_workflow.assert_called_once_with(other_workflow)

    def test_including_workflow_path_dispatches_to_include_path(self):
        workflow = Workflow()

        with patch.object(workflow, 'include_path') as mock_include_path:
            workflow.include('/path/to/other_workflow.py')
            mock_include_path.assert_called_once_with(
                '/path/to/other_workflow.py')

    @patch('gwf.core.inspect.ismodule', return_value=True)
    def test_including_workflow_module_gets_workflow_attribute_and_dispatches_to_include_workflow(self, mock_ismodule):
        workflow = Workflow(working_dir='/some/dir')
        other_workflow = Workflow(working_dir='/some/other/dir')

        mock_module = Mock()
        mock_module.gwf = other_workflow

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include(mock_module)

            mock_ismodule.assert_called_once_with(mock_module)
            mock_include_workflow.assert_called_once_with(other_workflow)

    def test_including_non_module_str_and_object_value_raises_type_error(self):
        workflow = Workflow(working_dir='/some/dir')
        with self.assertRaises(TypeError):
            workflow.include(42)

    def test_targets_inherit_workflow_working_dir_with_given_working_dir(self):
        workflow = Workflow(working_dir='/some/path')
        target = workflow.target('TestTarget', inputs=[], outputs=[])
        self.assertEqual(target.working_dir, '/some/path')

    @patch('gwf.core.sys._getframe')
    @patch('gwf.core.inspect.getfile', return_value='/some/path/file.py')
    def test_workflow_computes_working_dir_when_not_initialized_with_working_dir(
            self, inspect_getfile_mock, sys_getframe_mock):

        workflow = Workflow()

        self.assertEqual(sys_getframe_mock.call_count, 1)
        self.assertEqual(inspect_getfile_mock.call_count, 1)
        self.assertEqual(workflow.working_dir, '/some/path')

    def test_empty_workflow_repr(self):
        workflow = Workflow('/some/dir')
        self.assertEqual(
            repr(workflow), "Workflow(working_dir='/some/dir', targets={})")

    def test_nonempty_workflow_repr(self):
        workflow = Workflow('/some/dir')
        workflow.target('TestTarget')

        self.assertIn("working_dir='/some/dir'", repr(workflow))
        self.assertIn("name='TestTarget'", repr(workflow))


class TestTarget(unittest.TestCase):

    def test_relative_input_paths_are_normalized(self):
        target = create_test_target(
            inputs=['test_input1.txt', 'test_input2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(os.path.isabs(target.inputs[0]))
        self.assertTrue(os.path.isabs(target.inputs[1]))

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/some/path'))

    def test_relative_output_paths_are_normalized(self):
        target = create_test_target(
            outputs=['test_output1.txt', 'test_output2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(os.path.isabs(target.outputs[0]))
        self.assertTrue(os.path.isabs(target.outputs[1]))

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/some/path'))

    def test_absolute_input_paths_are_not_normalized(self):
        target = create_test_target(
            inputs=['test_input1.txt', '/other/path/test_input2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_absolute_output_paths_are_not_normalized(self):
        target = create_test_target(
            inputs=['test_output1.txt', '/other/path/test_output2.txt'],
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_target_without_outputs_is_a_sink(self):
        target = create_test_target()
        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = create_test_target(
            outputs=['test_output1.txt', 'test_output2.txt']
        )
        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = create_test_target()
        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = create_test_target(
            inputs=['test_input1.txt', 'test_input2.txt']
        )
        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = create_test_target() << 'this is a spec'
        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, 'this is a spec')

    def test_str_on_target(self):
        target = Target(
            'TestTarget',
            inputs=[],
            outputs=[],
            options={},
            working_dir='/some/dir'
        )
        self.assertEqual(str(target), 'TestTarget')

    def test_repr_on_target(self):
        target = Target(
            'TestTarget',
            inputs=[],
            outputs=[],
            options={},
            working_dir='/some/dir'
        )

        self.assertEqual(
            repr(target),
            "Target(name='TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir', spec=None)"
        )

    def test_repr_on_target_with_input(self):
        target = Target(
            'TestTarget',
            inputs=['test_input.txt'],
            outputs=[],
            options={},
            working_dir='/some/dir'
        )

        self.assertEqual(
            repr(target),
            "Target(name='TestTarget', inputs=['/some/dir/test_input.txt'], outputs=[], options={}, working_dir='/some/dir', spec=None)"
        )

    def test_repr_on_target_with_output(self):
        target = Target(
            'TestTarget',
            inputs=[],
            outputs=['test_output.txt'],
            options={},
            working_dir='/some/dir'
        )

        self.assertEqual(
            repr(target),
            "Target(name='TestTarget', inputs=[], outputs=['/some/dir/test_output.txt'], options={}, working_dir='/some/dir', spec=None)"
        )

    def test_repr_on_target_with_option(self):
        target = Target(
            'TestTarget',
            inputs=[],
            outputs=[],
            options={'some_option': True},
            working_dir='/some/dir'
        )

        self.assertEqual(
            repr(target),
            "Target(name='TestTarget', inputs=[], outputs=[], options={'some_option': True}, working_dir='/some/dir', spec=None)"
        )


class TestPreparedWorkflow(unittest.TestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')

    def test_finds_no_providers_in_empty_workflow(self):
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertDictEqual(prepared_workflow.provides, {})

    def test_finds_no_providers_in_workflow_with_no_producers(self):
        self.workflow.target('TestTarget', inputs=[], outputs=[])

        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertDictEqual(prepared_workflow.provides, {})

    def test_finds_provider_in_workflow_with_one_producer(self):
        self.workflow.target(
            'TestTarget', inputs=[], outputs=['/test_output.txt'], working_dir='')

        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertIn('/test_output.txt', prepared_workflow.provides)
        self.assertEqual(
            prepared_workflow.provides['/test_output.txt'].name, 'TestTarget')

    def test_raises_exceptions_if_two_targets_produce_the_same_file(self):
        self.workflow.target(
            'TestTarget1', inputs=[], outputs=['/test_output.txt'], working_dir='')
        self.workflow.target(
            'TestTarget2', inputs=[], outputs=['/test_output.txt'], working_dir='')

        with self.assertRaises(FileProvidedByMultipleTargetsError):
            PreparedWorkflow(self.workflow)

    def test_finds_no_dependencies_for_target_with_no_inputs(self):
        target = self.workflow.target('TestTarget', inputs=[], outputs=[])
        prepared_workflow = PreparedWorkflow(self.workflow)

        self.assertEqual(prepared_workflow.dependencies[target], [])

    @patch('gwf.core.os.path.exists', return_value=True)
    def test_existing_files_cause_no_dependencies(self, mock_os_path_exists):
        target = self.workflow.target(
            'TestTarget', inputs=['test_input.txt'], outputs=[])

        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertEqual(prepared_workflow.dependencies[target], [])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_non_existing_files_not_provided_by_other_target_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget', inputs=['test_input.txt'], outputs=[])
        with self.assertRaises(FileRequiredButNotProvidedError):
            PreparedWorkflow(self.workflow)

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_file_provided_by_other_target(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=['test_file.txt'], outputs=[])

        prepared_workflow = PreparedWorkflow(self.workflow)

        self.assertIn(target2, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target2])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_files_provided_by_two_other_targets(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file1.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=[], outputs=['test_file2.txt'])
        target3 = self.workflow.target(
            'TestTarget3', inputs=['test_file1.txt', 'test_file2.txt'], outputs=[])

        prepared_workflow = PreparedWorkflow(self.workflow)

        self.assertIn(target3, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target3])
        self.assertIn(target2, prepared_workflow.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_files_provided_by_other_targets_in_chain(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file1.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])
        target3 = self.workflow.target(
            'TestTarget3', inputs=['test_file2.txt'], outputs=[])

        prepared_workflow = PreparedWorkflow(self.workflow)

        self.assertIn(target2, prepared_workflow.dependencies)
        self.assertIn(target3, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target2])
        self.assertIn(target2, prepared_workflow.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_immediate_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget1', inputs=['test_file2.txt'], outputs=['test_file1.txt'])
        self.workflow.target(
            'TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])

        with self.assertRaises(CircularDependencyError):
            PreparedWorkflow(self.workflow)

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget1', inputs=['test_file3.txt'], outputs=['test_file1.txt'])
        self.workflow.target(
            'TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])
        self.workflow.target(
            'TestTarget3', inputs=['test_file2.txt'], outputs=['test_file3.txt'])

        with self.assertRaises(CircularDependencyError):
            PreparedWorkflow(self.workflow)

    def test_empty_prepared_workflow_repr(self):
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertEqual(
            repr(prepared_workflow),
            "PreparedWorkflow(working_dir='/some/dir', targets={})"
        )

    def test_nonempty_prepared_workflow_repr(self):
        self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.assertIn("working_dir='/some/dir'", repr(prepared_workflow))
        self.assertIn('TestTarget', repr(prepared_workflow))


class TestGetExecutionSchedule(unittest.TestCase):
    pass
