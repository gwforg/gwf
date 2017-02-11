import os
import os.path
import unittest
from unittest.mock import Mock, patch

import pathlib

from gwf import PreparedWorkflow, Target, Workflow, template, inputs, outputs, sink
from gwf.backends.testing import TestingBackend
from gwf.core import schedule, schedule_many
from gwf.exceptions import (CircularDependencyError,
                            FileProvidedByMultipleTargetsError,
                            FileRequiredButNotProvidedError,
                            IncludeWorkflowError, InvalidNameError,
                            TargetExistsError)

from . import create_test_target


class TestTemplate(unittest.TestCase):

    def setUp(self):
        self.test_template = template(
            inputs=['input{idx}.txt'],
            outputs=['output{idx}.txt'],
        )

        self.test_template << """
        cat input{idx}.txt > output{idx}.txt
        """

    def test_template_substitutes_args_when_called(self):

        options, spec = self.test_template(idx=0)
        self.assertDictEqual(options, {
            'inputs': ['input0.txt'],
            'outputs': ['output0.txt']
        })

        self.assertIn('cat input0.txt > output0.txt', spec)

    def test_target_created_from_template(self):
        workflow = Workflow(working_dir='/some/path')
        target = workflow.target('TestTarget') << self.test_template(idx=0)

        self.assertEqual(target.inputs, ['/some/path/input0.txt'])
        self.assertEqual(target.outputs, ['/some/path/output0.txt'])
        self.assertIn('cat input0.txt > output0.txt', target.spec)


class TestWorkflow(unittest.TestCase):

    def test_workflow_with_invalid_name_raises_error(self):
        with self.assertRaises(InvalidNameError):
            Workflow(name='123abc')

    def test_target_with_no_input_has_empty_inputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget') << sink()
        self.assertListEqual(target.inputs, [])

    def test_target_with_no_output_has_empty_outputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget') << sink()
        self.assertListEqual(target.outputs, [])

    def test_adding_a_target_makes_it_available_to_the_workflow(self):
        workflow = Workflow()
        workflow.target('TestTarget') << sink()

        self.assertIn('TestTarget', workflow.targets)

    def test_adding_two_targets_with_the_same_names_should_raise_an_exception(self):
        workflow = Workflow()
        workflow.target('TestTarget') << sink()

        with self.assertRaises(TargetExistsError):
            workflow.target('TestTarget') << sink()

    def test_including_workflow_with_no_name_raises_an_exception(self):
        workflow = Workflow()
        other_workflow = Workflow()
        with self.assertRaises(IncludeWorkflowError):
            workflow.include(other_workflow)

    def test_including_workflow_object_should_extend_including_workflow(self):
        workflow = Workflow()
        workflow.target('TestTarget1') << sink()

        other_workflow = Workflow(name='foo')
        other_workflow.target('TestTarget2') << sink()
        other_workflow.target('TestTarget3') << sink()

        workflow.include_workflow(other_workflow)

        self.assertIn('TestTarget1', workflow.targets)
        self.assertIn('foo.TestTarget2', workflow.targets)
        self.assertIn('foo.TestTarget3', workflow.targets)

    def test_include_with_namespace_overrides_included_workflow_name(self):
        workflow = Workflow()
        workflow.target('TestTarget') << sink()

        other_workflow = Workflow(name='foo')
        other_workflow.target('TestTarget') << sink()

        workflow.include_workflow(other_workflow, namespace='bar')
        self.assertIn('bar.TestTarget', workflow.targets)
        self.assertNotIn('foo.TestTarget', workflow.targets)
        self.assertIn('TestTarget', workflow.targets)

    def test_including_workflow_with_same_name_as_this_workflow_raises_an_exception(self):
        workflow = Workflow(name='foo')
        other_workflow = Workflow(name='foo')
        with self.assertRaises(IncludeWorkflowError):
            workflow.include(other_workflow)

    @patch('gwf.core.import_object')
    def test_including_workflow_path_import_object_and_include_workflow_into_current_workflow(self, mock_import_object):
        workflow = Workflow()
        workflow.target('TestTarget1') << sink()

        other_workflow = Workflow()
        other_workflow.target('TestTarget2') << sink()
        other_workflow.target('TestTarget3') << sink()

        mock_import_object.return_value = other_workflow

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include_path('/path/to/other_workflow.py')

            mock_import_object.assert_called_once_with(
                '/path/to/other_workflow.py'
            )

            mock_include_workflow.assert_called_once_with(
                other_workflow, namespace=None
            )

    def test_including_workflow_instance_dispatches_to_include_workflow(self):
        workflow = Workflow()
        other_workflow = Workflow()

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include(other_workflow)
            mock_include_workflow.assert_called_once_with(
                other_workflow, namespace=None)

    def test_including_workflow_path_dispatches_to_include_path(self):
        workflow = Workflow()

        with patch.object(workflow, 'include_path') as mock_include_path:
            workflow.include('/path/to/other_workflow.py')
            mock_include_path.assert_called_once_with(
                '/path/to/other_workflow.py', namespace=None)

    @patch('gwf.core.inspect.ismodule', return_value=True)
    def test_including_workflow_module_gets_workflow_attribute_and_dispatches_to_include_workflow(self, mock_ismodule):
        workflow = Workflow(working_dir='/some/dir')
        other_workflow = Workflow(working_dir='/some/other/dir')

        mock_module = Mock()
        mock_module.gwf = other_workflow

        with patch.object(workflow, 'include_workflow') as mock_include_workflow:
            workflow.include(mock_module)

            mock_ismodule.assert_called_once_with(mock_module)
            mock_include_workflow.assert_called_once_with(
                other_workflow, namespace=None)

    def test_including_non_module_str_and_object_value_raises_type_error(self):
        workflow = Workflow(working_dir='/some/dir')
        with self.assertRaises(TypeError):
            workflow.include(42)

    def test_targets_inherit_workflow_working_dir_with_given_working_dir(self):
        workflow = Workflow(working_dir='/some/path')
        target = workflow.target('TestTarget') << sink()
        self.assertEqual(target.working_dir, '/some/path')

    @patch('gwf.core.sys._getframe')
    @patch('gwf.core.inspect.getfile', return_value='/some/path/file.py')
    def test_workflow_computes_working_dir_when_not_initialized_with_working_dir(
            self, inspect_getfile_mock, sys_getframe_mock):

        workflow = Workflow()

        self.assertEqual(sys_getframe_mock.call_count, 1)
        self.assertEqual(inspect_getfile_mock.call_count, 1)
        self.assertEqual(workflow.working_dir, '/some/path')

    @patch('gwf.core._glob')
    def test_glob_with_relative_path_searches_relative_to_working_dir(self, glob_mock):
        workflow = Workflow(working_dir='/some/path')
        workflow.glob('*.fa')
        glob_mock.assert_called_once_with('/some/path/*.fa')

    @patch('gwf.core._glob', return_value=['/other/path/A.fa', '/other/path/B.fa'])
    def test_glob_with_absolute_path_does_not_search_relative_to_working_dir(self, glob_mock):
        workflow = Workflow(working_dir='/some/path')
        res = workflow.glob('/other/path/*.fa')
        self.assertEqual(res, ['/other/path/A.fa', '/other/path/B.fa'])
        glob_mock.assert_called_once_with('/other/path/*.fa')

    @patch('gwf.core._iglob')
    def test_iglob_with_relative_path_searches_relative_to_working_dir(self, iglob_mock):
        workflow = Workflow(working_dir='/some/path')
        workflow.iglob('*.fa')
        iglob_mock.assert_called_once_with('/some/path/*.fa')

    @patch('gwf.core._iglob', return_value=['/other/path/A.fa', '/other/path/B.fa'])
    def test_iglob_with_absolute_path_does_not_search_relative_to_working_dir(self, iglob_mock):
        workflow = Workflow(working_dir='/some/path')
        res = list(workflow.iglob('/other/path/*.fa'))
        self.assertEqual(res, ['/other/path/A.fa', '/other/path/B.fa'])
        iglob_mock.assert_called_once_with('/other/path/*.fa')

    @patch('gwf.core.subprocess.check_output')
    def test_shell_calls_subprocess_with_same_working_dir_as_workflow_in_a_shell(self, mock_check_output):
        workflow = Workflow(working_dir='/some/path')
        workflow.shell('echo hello')
        mock_check_output.assert_called_once_with('echo hello', cwd='/some/path', shell=True)


class TestTarget(unittest.TestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/path')

    def test_target_with_invalid_name_raises_exception(self):
        with self.assertRaises(InvalidNameError):
            Target('123abc', workflow=self.workflow) << sink()

    def test_relative_input_paths_are_normalized(self):
        target = create_test_target(
            input_files=['test_input1.txt', 'test_input2.txt'],
            workflow=self.workflow
        )

        self.assertTrue(os.path.isabs(target.inputs[0]))
        self.assertTrue(os.path.isabs(target.inputs[1]))

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/some/path'))

    def test_relative_output_paths_are_normalized(self):
        target = create_test_target(
            output_files=['test_output1.txt', 'test_output2.txt'],
            workflow=self.workflow
        )

        self.assertTrue(os.path.isabs(target.outputs[0]))
        self.assertTrue(os.path.isabs(target.outputs[1]))

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/some/path'))

    def test_absolute_input_paths_are_not_normalized(self):
        target = create_test_target(
            input_files=['test_input1.txt', '/other/path/test_input2.txt'],
            workflow=self.workflow
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_absolute_output_paths_are_not_normalized(self):
        target = create_test_target(
            output_files=['test_output1.txt', '/other/path/test_output2.txt'],
            workflow=self.workflow
        )

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/other/path'))

    def test_target_without_outputs_is_a_sink(self):
        target = create_test_target()
        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = create_test_target(
            output_files=['test_output1.txt', 'test_output2.txt']
        )
        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = create_test_target()
        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = create_test_target(
            input_files=['test_input1.txt', 'test_input2.txt']
        )
        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = create_test_target() << 'this is a spec'
        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, 'this is a spec')

    def test_raises_valueerror_if_inputs_is_not_valid(self):
        pass  # This is no longer an error...
        # with self.assertRaises(InvalidTypeError):
        #     create_test_target(input_files='hello.txt')

    def test_raises_valueerror_if_outputs_is_not_valid(self):
        pass  # This is no longer an error...
        # with self.assertRaises(InvalidTypeError):
        #     create_test_target(output_files='hello.txt')

    def test_should_stringify_input_paths(self):
        target = create_test_target(input_files=[pathlib.PurePath('somefile.txt'), 'otherfile.txt'])
        self.assertEqual(target.inputs, ['/some/path/somefile.txt', '/some/path/otherfile.txt'])

    def test_should_stringify_output_paths(self):
        target = create_test_target(output_files=[pathlib.PurePath('somefile.txt'), 'otherfile.txt'])
        self.assertEqual(target.outputs, ['/some/path/somefile.txt', '/some/path/otherfile.txt'])

    def test_warn_user_that_target_is_a_sink_when_warn_sink_is_true(self):
        with self.assertLogs(level='WARNING') as logs:
            target = create_test_target(warn_sink=True)
            PreparedWorkflow(
                targets={'target': target},
                working_dir='/some/dir',
                supported_options={},
                config={},
            )

            self.assertEqual(len(logs.output), 1)
            self.assertIn('TestTarget', logs.output[0])
            self.assertIn('sink()', logs.output[0])

    def test_do_not_warn_user_that_target_is_a_sink_when_warn_sink_is_false(self):
        try:
            # This is weird. When accessing `logs.output` and nothing was logged,
            # an AssertionError is raised saying nothing was logged. Thus, we need
            # to catch this error. If no error was raised, something must have
            # logged and thus we fail the test.
            with self.assertLogs(level='WARNING') as logs:
                create_test_target(warn_sink=False)
                self.assertEqual(logs.output, [])
        except AssertionError:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_str_on_target(self):
        target = Target('TestTarget', workflow=self.workflow)
        self.assertEqual(str(target), 'TestTarget')


class TestPreparedWorkflow(unittest.TestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.supported_options = {}
        self.config = {}

    def test_finds_no_providers_in_empty_workflow(self):
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )
        self.assertDictEqual(prepared_workflow.provides, {})

    def test_finds_no_providers_in_workflow_with_no_producers(self):
        self.workflow.target('TestTarget') << sink()

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )
        self.assertDictEqual(prepared_workflow.provides, {})

    def test_finds_provider_in_workflow_with_one_producer(self):
        self.workflow.target('TestTarget') <<\
            outputs('/test_output.txt')

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )
        self.assertIn('/test_output.txt', prepared_workflow.provides)
        self.assertEqual(
            prepared_workflow.provides['/test_output.txt'].name, 'TestTarget')

    def test_raises_exceptions_if_two_targets_produce_the_same_file(self):
        self.workflow.target('TestTarget1') <<\
            outputs('/test_output.txt')
        self.workflow.target('TestTarget2') <<\
            outputs('/test_output.txt')

        with self.assertRaises(FileProvidedByMultipleTargetsError):
            PreparedWorkflow(
                targets=self.workflow.targets,
                working_dir=self.workflow.working_dir,
                supported_options=self.supported_options,
                config=self.config,
            )

    def test_finds_no_dependencies_for_target_with_no_inputs(self):
        target = self.workflow.target('TestTarget') << sink()
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )

        self.assertEqual(prepared_workflow.dependencies[target], [])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_non_existing_files_not_provided_by_other_target_raises_exception(self, mock_os_path_exists):
        self.workflow.target('TestTarget') << inputs('test_input.txt')
        with self.assertRaises(FileRequiredButNotProvidedError):
            PreparedWorkflow(
                targets=self.workflow.targets,
                working_dir=self.workflow.working_dir,
                supported_options=self.supported_options,
                config=self.config,
            )

    @patch('gwf.core.os.path.exists', return_value=True)
    def test_existing_files_not_provided_by_other_target_has_no_dependencies(self, mock_exists):
        target = self.workflow.target('TestTarget') << inputs('test_file.txt') << sink()
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )
        self.assertListEqual(prepared_workflow.dependencies[target], [])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_file_provided_by_other_target(self, mock_os_path_exists):
        target1 = self.workflow.target('TestTarget1') <<\
            outputs('test_file.txt')
        target2 = self.workflow.target('TestTarget2') <<\
            inputs('test_file.txt') << sink()

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )

        self.assertIn(target2, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target2])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_files_provided_by_two_other_targets(self, mock_os_path_exists):
        target1 = self.workflow.target('TestTarget1') <<\
            outputs('test_file1.txt')
        target2 = self.workflow.target('TestTarget2') <<\
            outputs('test_file2.txt')
        target3 = self.workflow.target('TestTarget3') <<\
            inputs('test_file1.txt', 'test_file2.txt') <<\
            sink()

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )

        self.assertIn(target3, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target3])
        self.assertIn(target2, prepared_workflow.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_finds_non_existing_files_provided_by_other_targets_in_chain(self, mock_os_path_exists):
        target1 = self.workflow.target('TestTarget1') <<\
            outputs('test_file1.txt')
        target2 = self.workflow.target('TestTarget2') <<\
            inputs('test_file1.txt') <<\
            outputs('test_file2.txt')
        target3 = self.workflow.target('TestTarget3') <<\
            inputs('test_file2.txt') << sink()

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )

        self.assertIn(target2, prepared_workflow.dependencies)
        self.assertIn(target3, prepared_workflow.dependencies)
        self.assertIn(target1, prepared_workflow.dependencies[target2])
        self.assertIn(target2, prepared_workflow.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_immediate_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget1') <<\
            inputs('test_file2.txt') <<\
            outputs('test_file1.txt')
        self.workflow.target('TestTarget2') <<\
            inputs('test_file1.txt') <<\
            outputs('test_file2.txt')

        with self.assertRaises(CircularDependencyError):
            PreparedWorkflow(
                targets=self.workflow.targets,
                working_dir=self.workflow.working_dir,
                supported_options=self.supported_options,
                config=self.config,
            )

    @patch('gwf.core.os.path.exists', return_value=False)
    def test_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget1') <<\
            inputs('test_file3.txt') <<\
            outputs('test_file1.txt')
        self.workflow.target(
            'TestTarget2') <<\
            inputs('test_file1.txt') <<\
            outputs('test_file2.txt')
        self.workflow.target(
            'TestTarget3') <<\
            inputs('test_file2.txt') <<\
            outputs('test_file3.txt')

        with self.assertRaises(CircularDependencyError):
            PreparedWorkflow(
                targets=self.workflow.targets,
                working_dir=self.workflow.working_dir,
                supported_options=self.supported_options,
                config=self.config,
            )

    def test_endpoints_only_includes_target_with_no_dependents(self):
        self.workflow.target('TestTarget1') << outputs('test.txt')
        target2 = self.workflow.target('TestTarget2') << inputs('test.txt') << sink()
        target3 = self.workflow.target('TestTarget3') << inputs('test.txt') << sink()
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=self.supported_options,
            config=self.config,
        )
        self.assertSetEqual(prepared_workflow.endpoints(), {target2, target3})


class TestShouldRun(unittest.TestCase):

    def setUp(self):
        workflow = Workflow(working_dir='/some/dir')
        self.target1 = workflow.target('TestTarget1') <<\
            outputs('test_output1.txt')
        self.target2 = workflow.target('TestTarget2') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output2.txt')
        self.target3 = workflow.target('TestTarget3') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output3.txt')
        self.target4 = workflow.target('TestTarget4') <<\
            inputs('test_output2.txt', 'test_output3.txt') <<\
            outputs('final_output.txt')

        self.prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options={},
            config={},
        )

    def test_target_should_run_if_one_of_its_dependencies_does_not_exist(self):
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(self.prepared_workflow.should_run(self.target1))

        self.assertEqual(
            logs.output,
            [
                'DEBUG:gwf.core:TestTarget1 should run because one of its output files does not exist.'
            ]
        )

    def test_target_should_run_if_one_of_its_dependencies_should_run(self):
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(self.prepared_workflow.should_run(self.target2))

        self.assertEqual(
            logs.output,
            [
                'DEBUG:gwf.core:TestTarget1 should run because one of its output files does not exist.',
                'DEBUG:gwf.core:TestTarget2 should run because one of its dependencies should run.'
            ]
        )

    def test_target_should_run_if_any_of_its_deep_dependencies_should_run(self):
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(self.prepared_workflow.should_run(self.target4))

        self.assertEqual(
            logs.output,
            [
                'DEBUG:gwf.core:TestTarget1 should run because one of its output files does not exist.',
                'DEBUG:gwf.core:TestTarget2 should run because one of its dependencies should run.',
                'DEBUG:gwf.core:TestTarget4 should run because one of its dependencies should run.'
            ]
        )

    def test_target_should_run_if_it_is_a_sink(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget') << sink()

        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options={},
            config={},
        )

        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(prepared_workflow.should_run(target))

            self.assertEqual(
                logs.output,
                [
                    'DEBUG:gwf.core:TestTarget should run because it is a sink.'
                ]
            )

    def test_target_should_not_run_if_it_is_a_source_and_all_outputs_exist(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget1') <<\
            outputs('test_output1.txt', 'test_output2.txt')

        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options={},
            config={},
        )

        mock_file_cache = {
            '/some/dir/test_output1.txt': 1,
            '/some/dir/test_output2.txt': 2
        }
        with patch.dict(prepared_workflow.file_cache, mock_file_cache):
            self.assertFalse(
                prepared_workflow.should_run(target)
            )

    def test_should_run_if_any_input_file_is_newer_than_any_output_file(self):
        mock_file_cache = {
            '/some/dir/test_output1.txt': 0,
            '/some/dir/test_output2.txt': 1,
            '/some/dir/test_output3.txt': 3,
            '/some/dir/final_output.txt': 2,
        }

        with patch.dict(self.prepared_workflow.file_cache, mock_file_cache):
            self.assertFalse(self.prepared_workflow.should_run(self.target1))
            self.assertFalse(self.prepared_workflow.should_run(self.target2))
            self.assertFalse(self.prepared_workflow.should_run(self.target3))
            self.assertTrue(self.prepared_workflow.should_run(self.target4))

    def test_should_run_not_run_if_all_outputs_are_newer_then_the_inputs(self):
        mock_file_cache = {
            '/some/dir/test_output1.txt': 0,
            '/some/dir/test_output2.txt': 1,
            '/some/dir/test_output3.txt': 3,
            '/some/dir/final_output.txt': 4,
        }

        with patch.dict(self.prepared_workflow.file_cache, mock_file_cache):
            self.assertFalse(self.prepared_workflow.should_run(self.target1))
            self.assertFalse(self.prepared_workflow.should_run(self.target2))
            self.assertFalse(self.prepared_workflow.should_run(self.target3))
            self.assertFalse(self.prepared_workflow.should_run(self.target4))

    @patch('gwf.core.os.path.exists', return_value=True)
    def test_two_targets_producing_the_same_file_but_declared_with_rel_and_abs_path(self, mock_os_path_exists):
        workflow = Workflow(working_dir='/some/dir')
        workflow.target('TestTarget1') << outputs('/some/dir/test_output.txt')
        workflow.target('TestTarget2') << outputs('test_output.txt')

        with self.assertRaises(FileProvidedByMultipleTargetsError):
            PreparedWorkflow(
                targets=workflow.targets,
                working_dir=workflow.working_dir,
                supported_options={},
                config={},
            )


class TestScheduling(unittest.TestCase):

    def test_scheduling_workflow_with_one_target_that_is_already_submitted_returns_empty_schedule(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget') << sink()

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=True):
            sched = schedule(prepared_workflow, backend, target)
            self.assertListEqual(sched, [])

    def test_scheduling_workflow_with_one_target_that_is_not_submitted_returns_schedule_with_target(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget') << sink()

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                sched = schedule(prepared_workflow, backend, target)
                self.assertListEqual(sched, [target])

    def test_scheduling_workflow_with_target_with_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1') << outputs('test_output.txt')
        target2 = workflow.target('TestTarget2') <<\
            inputs('test_output.txt') <<\
            sink()

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                sched = schedule(prepared_workflow, backend, target2)
                self.assertListEqual(sched, [target1, target2])

    def test_scheduling_workflow_with_deep_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1') << outputs('test_output1.txt')
        target2 = workflow.target('TestTarget2') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output2.txt')
        target3 = workflow.target('TestTarget3') <<\
            inputs('test_output2.txt') <<\
            outputs('test_output3.txt')
        target4 = workflow.target('TestTarget4') <<\
            inputs('test_output3.txt') <<\
            outputs('final_output.txt')

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                sched = schedule(prepared_workflow, backend, target4)
                self.assertListEqual(
                    sched, [target1, target2, target3, target4])

    def test_scheduling_workflow_with_branch_and_join_structure(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1') << outputs('test_output1.txt')
        target2 = workflow.target('TestTarget2') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output2.txt')
        target3 = workflow.target('TestTarget3') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output3.txt')
        target4 = workflow.target('TestTarget4') <<\
            inputs('test_output2.txt', 'test_output3.txt') <<\
            outputs('final_output.txt')

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                sched = schedule(prepared_workflow, backend, target4)
                self.assertListEqual(
                    sched, [target1, target2, target3, target4])

    def test_option_defaults_contains_specified_option_defaults(self):
        backend = TestingBackend()
        backend.configure(working_dir='/some/dir', config={})
        self.assertDictEqual(
            backend.option_defaults,
            {'cores': 2, 'memory': '18gb'}
        )

    def test_scheduling_a_submitted_dependency_does_not_submit_the_dependency(self):
        workflow = Workflow(working_dir='/some/dir')

        workflow.target('TestTarget1') <<\
            outputs('test_output1.txt')
        target2 = workflow.target('TestTarget2') <<\
            outputs('test_output2.txt')
        target3 = workflow.target('TestTarget3') <<\
            inputs('test_output1.txt', 'test_output2.txt') <<\
            outputs('test_output3.txt')

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', side_effect=[False, True, False, False]):
            sched = schedule(prepared_workflow, backend, target3)
            self.assertEqual(sched, [target2, target3])

    def test_scheduling_non_submitted_targets_that_should_not_run_does_not_submit_any_targets(self):
        workflow = Workflow(working_dir='/some/dir')

        workflow.target('TestTarget1') << outputs('test_output1.txt')
        workflow.target('TestTarget2') << outputs('test_output2.txt')
        target3 = workflow.target('TestTarget3') <<\
            inputs('test_output1.txt', 'test_output2.txt') <<\
            outputs('test_output3.txt')

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', side_effect=[False, False, False, False]):
                sched = schedule(prepared_workflow, backend, target3)
                self.assertEqual(sched, [])

    def test_scheduling_many_targets_calls_schedule_for_each_target(self):
        workflow = Workflow(working_dir='/some/dir')

        target1 = workflow.target('TestTarget1') << outputs('test_output1.txt')
        target2 = workflow.target('TestTarget2') << outputs('test_output2.txt')
        target3 = workflow.target('TestTarget3') <<\
            inputs('test_output1.txt') <<\
            outputs('test_output3.txt')
        target4 = workflow.target('TestTarget4') <<\
            inputs('test_output2.txt') <<\
            outputs('test_output4.txt')

        backend = TestingBackend()
        prepared_workflow = PreparedWorkflow(
            targets=workflow.targets,
            working_dir=workflow.working_dir,
            supported_options=backend.supported_options,
            config={},
        )
        backend.configure(working_dir=prepared_workflow.working_dir, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                sched = schedule_many(
                    prepared_workflow, backend, [target3, target4])
                self.assertEqual(
                    sched, [[target1, target3], [target2, target4]])
