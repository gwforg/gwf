import os
import os.path
import pathlib
import unittest
from unittest.mock import Mock, patch, call

import pytest

from gwf import Graph, Target, Scheduler, Workflow
from gwf.backends import Backend, Status
from gwf.exceptions import (CircularDependencyError,
                            MultipleProvidersError,
                            MissingProviderError,
                            IncludeWorkflowError, InvalidNameError,
                            TargetExistsError, InvalidTypeError)
from gwf.backends.exceptions import LogNotFoundError


class DummyBackend(Backend):
    def __init__(self):
        super().__init__()
        self._tracked = {}

    def submit(self, target, dependencies):
        self._tracked[target] = Status.SUBMITTED

    def cancel(self, target):
        del self._tracked[target]

    def status(self, target):
        return self._tracked.get(target, Status.UNKNOWN)

    def logs(self, target, stderr=False):
        raise LogNotFoundError

    def close(self):
        pass

    def set_status(self, target, status):
        assert status in (Status.RUNNING, Status.UNKNOWN)
        self._tracked[target] = status


@pytest.fixture
def backend():
    backend = DummyBackend()
    backend.submit = Mock(wraps=backend.submit)
    backend.status = Mock(wraps=backend.status)
    return backend


class TestWorkflow(unittest.TestCase):
    def test_workflow_with_invalid_name_raises_error(self):
        with self.assertRaises(InvalidNameError):
            Workflow(name='123abc')

    def test_target_with_no_input_has_empty_inputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget', inputs=[], outputs=[])
        self.assertListEqual(target.inputs, [])

    def test_target_with_no_output_has_empty_outputs_attribute(self):
        workflow = Workflow()
        target = workflow.target('TestTarget', inputs=[], outputs=[])
        self.assertListEqual(target.outputs, [])

    def test_adding_a_target_makes_it_available_to_the_workflow(self):
        workflow = Workflow()
        target = workflow.target('TestTarget', inputs=[], outputs=[])

        self.assertIn('TestTarget', workflow.targets)
        self.assertIn(target, workflow.targets.values())

    def test_adding_two_targets_with_the_same_names_should_raise_an_exception(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=[])

        with self.assertRaises(TargetExistsError):
            workflow.target('TestTarget', inputs=[], outputs=[])

    def test_including_workflow_with_no_name_raises_an_exception(self):
        workflow = Workflow()
        other_workflow = Workflow()
        with self.assertRaises(IncludeWorkflowError):
            workflow.include(other_workflow)

    def test_including_workflow_object_should_extend_including_workflow(self):
        workflow = Workflow()
        workflow.target('TestTarget1', inputs=[], outputs=[])

        other_workflow = Workflow(name='foo')
        other_workflow.target('TestTarget2', inputs=[], outputs=[])
        other_workflow.target('TestTarget3', inputs=[], outputs=[])

        workflow.include_workflow(other_workflow)

        self.assertIn('TestTarget1', workflow.targets)
        self.assertIn('foo.TestTarget2', workflow.targets)
        self.assertIn('foo.TestTarget3', workflow.targets)

    def test_include_target_from_workflow_in_two_different_workflows_(self):
        w1 = Workflow()
        target = w1.target('MyTarget', inputs=[], outputs=[])

        w3 = Workflow()
        w3.include(w1, namespace='bar')

        w2 = Workflow()
        w2.include(w1, namespace='foo')

        self.assertEqual(target.name, 'MyTarget')
        self.assertIn('bar.MyTarget', w3.targets)
        self.assertEqual(w3.targets['bar.MyTarget'].name, 'bar.MyTarget')
        self.assertIn('foo.MyTarget', w2.targets)
        self.assertEqual(w2.targets['foo.MyTarget'].name, 'foo.MyTarget')

    def test_including_workflow_with_same_name_as_this_workflow_raises_an_exception(self):
        workflow = Workflow(name='foo')
        other_workflow = Workflow(name='foo')
        with self.assertRaises(IncludeWorkflowError):
            workflow.include(other_workflow)

    @patch('gwf.core.load_workflow', autospec=True)
    def test_including_workflow_from_path(self, mock_load_workflow):
        workflow = Workflow()
        target1 = workflow.target('TestTarget1', inputs=[], outputs=[])

        other_workflow = Workflow()
        target2 = other_workflow.target('TestTarget2', inputs=[], outputs=[])
        target3 = other_workflow.target('TestTarget3', inputs=[], outputs=[])

        mock_load_workflow.return_value = other_workflow

        workflow.include_path('/path/to/other_workflow.py', namespace='other')
        self.assertEqual(workflow.targets.keys(), {'TestTarget1', 'other.TestTarget2', 'other.TestTarget3'})

    def test_including_workflow_instance_dispatches_to_include_workflow(self):
        workflow = Workflow()
        other_workflow = Workflow()

        with patch.object(workflow, 'include_workflow', autospec=True) as mock_include_workflow:
            workflow.include(other_workflow)
            mock_include_workflow.assert_called_once_with(
                other_workflow, namespace=None)

    def test_including_workflow_path_dispatches_to_include_path(self):
        workflow = Workflow()

        with patch.object(workflow, 'include_path', autospec=True) as mock_include_path:
            workflow.include('/path/to/other_workflow.py')
            mock_include_path.assert_called_once_with(
                '/path/to/other_workflow.py', namespace=None)

    @patch('gwf.core.inspect.ismodule', return_value=True, autospec=True)
    def test_including_workflow_module_gets_workflow_attribute_and_dispatches_to_include_workflow(self, mock_ismodule):
        workflow = Workflow(working_dir='/some/dir')
        other_workflow = Workflow(working_dir='/some/other/dir')

        mock_module = Mock()
        mock_module.gwf = other_workflow

        with patch.object(workflow, 'include_workflow', autospec=True) as mock_include_workflow:
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
        target = workflow.target('TestTarget', inputs=[], outputs=[])
        self.assertEqual(target.working_dir, '/some/path')

    def test_targets_inherit_workflow_defaults(self):
        workflow = Workflow(defaults={'cores': 8, 'memory': '8g'})
        target = workflow.target('TestTarget', inputs=[], outputs=[])
        self.assertEqual(target.options, {'cores': 8, 'memory': '8g'})

    def test_target_options_override_defaults(self):
        workflow = Workflow(defaults={'cores': 8, 'memory': '8g'})
        target = workflow.target('TestTarget', inputs=[], outputs=[], cores=16)
        self.assertEqual(target.options, {'cores': 16, 'memory': '8g'})

    @patch('gwf.core.sys._getframe', autospec=True)
    @patch('gwf.core.inspect.getfile', return_value='/some/path/file.py', autospec=True)
    def test_workflow_computes_working_dir_when_not_initialized_with_working_dir(
            self, inspect_getfile_mock, sys_getframe_mock):
        workflow = Workflow()

        self.assertEqual(sys_getframe_mock.call_count, 1)
        self.assertEqual(inspect_getfile_mock.call_count, 1)
        self.assertEqual(workflow.working_dir, '/some/path')

    @patch('gwf.core._glob', autospec=True)
    def test_glob_with_relative_path_searches_relative_to_working_dir(self, glob_mock):
        workflow = Workflow(working_dir='/some/path')
        workflow.glob('*.fa')
        glob_mock.assert_called_once_with('/some/path/*.fa')

    @patch('gwf.core._glob', return_value=['/other/path/A.fa', '/other/path/B.fa'], autospec=True)
    def test_glob_with_absolute_path_does_not_search_relative_to_working_dir(self, glob_mock):
        workflow = Workflow(working_dir='/some/path')
        res = workflow.glob('/other/path/*.fa')
        self.assertEqual(res, ['/other/path/A.fa', '/other/path/B.fa'])
        glob_mock.assert_called_once_with('/other/path/*.fa')

    @patch('gwf.core._iglob', autospec=True)
    def test_iglob_with_relative_path_searches_relative_to_working_dir(self, iglob_mock):
        workflow = Workflow(working_dir='/some/path')
        workflow.iglob('*.fa')
        iglob_mock.assert_called_once_with('/some/path/*.fa')

    @patch('gwf.core._iglob', return_value=['/other/path/A.fa', '/other/path/B.fa'], autospec=True)
    def test_iglob_with_absolute_path_does_not_search_relative_to_working_dir(self, iglob_mock):
        workflow = Workflow(working_dir='/some/path')
        res = list(workflow.iglob('/other/path/*.fa'))
        self.assertEqual(res, ['/other/path/A.fa', '/other/path/B.fa'])
        iglob_mock.assert_called_once_with('/other/path/*.fa')

    @patch('gwf.core.subprocess.check_output', autospec=True)
    def test_shell_calls_subprocess_with_same_working_dir_as_workflow_in_a_shell(self, mock_check_output):
        workflow = Workflow(working_dir='/some/path')
        workflow.shell('echo hello')
        mock_check_output.assert_called_once_with('echo hello', cwd='/some/path', shell=True)


class TestTarget(unittest.TestCase):
    def test_target_with_invalid_name_raises_exception(self):
        with self.assertRaises(InvalidNameError):
            Target(
                '123abc',
                inputs=[],
                outputs=[],
                options={},
                working_dir='/some/path'
            )

    def test_relative_input_paths_are_normalized(self):
        target = Target(
            name='TestTarget',
            inputs=['test_input1.txt', 'test_input2.txt'],
            outputs=[],
            options={},
            working_dir='/some/path',
        )

        self.assertTrue(os.path.isabs(target.inputs[0]))
        self.assertTrue(os.path.isabs(target.inputs[1]))

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/some/path'))

    def test_relative_output_paths_are_normalized(self):
        target = Target(
            name='TestTarget',
            inputs=[],
            outputs=['test_output1.txt', 'test_output2.txt'],
            options={},
            working_dir='/some/path'
        )

        self.assertTrue(os.path.isabs(target.outputs[0]))
        self.assertTrue(os.path.isabs(target.outputs[1]))

        self.assertTrue(target.outputs[0].startswith('/some/path'))
        self.assertTrue(target.outputs[1].startswith('/some/path'))

    def test_absolute_input_paths_are_not_normalized(self):
        target = Target(
            name='TestTarget',
            inputs=['test_input1.txt', '/other/path/test_input2.txt'],
            outputs=[],
            options={},
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_absolute_output_paths_are_not_normalized(self):
        target = Target(
            name='TestTarget',
            inputs=['test_output1.txt', '/other/path/test_output2.txt'],
            outputs=[],
            options={},
            working_dir='/some/path'
        )

        self.assertTrue(target.inputs[0].startswith('/some/path'))
        self.assertTrue(target.inputs[1].startswith('/other/path'))

    def test_target_without_outputs_is_a_sink(self):
        target = Target(name='TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/path')
        self.assertTrue(target.is_sink)

    def test_target_with_outputs_is_not_a_sink(self):
        target = Target(
            name='TestTarget',
            inputs=[],
            outputs=['test_output1.txt', 'test_output2.txt'],
            options={},
            working_dir='/some/path',
        )
        self.assertFalse(target.is_sink)

    def test_target_without_inputs_is_a_source(self):
        target = Target(name='TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/path')
        self.assertTrue(target.is_source)

    def test_target_with_inputs_is_not_a_source(self):
        target = Target(
            name='TestTarget',
            inputs=['test_input1.txt', 'test_input2.txt'],
            outputs=[],
            options={},
            working_dir='/some/path'
        )
        self.assertFalse(target.is_source)

    def test_assigning_spec_to_target_sets_spec_attribute(self):
        target = Target(name='TestTarget', inputs=[], outputs=[], options={},
                        working_dir='/some/path') << 'this is a spec'
        self.assertIsNotNone(target.spec)
        self.assertEqual(target.spec, 'this is a spec')

    def test_raises_valueerror_if_inputs_is_not_valid(self):
        with self.assertRaises(InvalidTypeError):
            Target(name='TestTarget', inputs='hello.txt', outputs=[], options={}, working_dir='/some/path')

    def test_raises_valueerror_if_outputs_is_not_valid(self):
        with self.assertRaises(InvalidTypeError):
            Target(name='TestTarget', inputs=[], outputs='hello.txt', options={}, working_dir='/some/path')

    def test_should_stringify_input_paths(self):
        target = Target(
            name='TestTarget',
            inputs=[pathlib.PurePath('somefile.txt'), 'otherfile.txt'],
            outputs=[],
            options={},
            working_dir='/some/path'
        )
        self.assertEqual(target.inputs, ['/some/path/somefile.txt', '/some/path/otherfile.txt'])

    def test_should_stringify_output_paths(self):
        target = Target(
            name='TestTarget',
            inputs=[],
            outputs=[pathlib.PurePath('somefile.txt'), 'otherfile.txt'],
            options={},
            working_dir='/some/path'
        )
        self.assertEqual(target.outputs, ['/some/path/somefile.txt', '/some/path/otherfile.txt'])

    def test_inherit_options(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={'cores': 8}, working_dir='/some/dir')
        target.inherit_options({'cores': 4, 'memory': '4g'})
        self.assertEqual(target.options, {'cores': 8, 'memory': '4g'})

    def test_str_on_target(self):
        target = Target(
            'TestTarget',
            inputs=[],
            outputs=[],
            options={},
            working_dir='/some/path'
        )
        self.assertEqual(str(target), 'TestTarget')


class TestGraph(unittest.TestCase):
    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.supported_options = {}
        self.config = {}

    def test_finds_no_providers_in_empty_workflow(self):
        graph = Graph.from_targets(self.workflow.targets)
        self.assertDictEqual(graph.provides, {})

    def test_finds_no_providers_in_workflow_with_no_producers(self):
        self.workflow.target('TestTarget', inputs=[], outputs=[])
        graph = Graph.from_targets(self.workflow.targets)
        self.assertDictEqual(graph.provides, {})

    def test_finds_provider_in_workflow_with_one_producer(self):
        self.workflow.target('TestTarget', inputs=[], outputs=['/test_output.txt'], working_dir='')
        graph = Graph.from_targets(self.workflow.targets)
        self.assertIn('/test_output.txt', graph.provides)
        self.assertEqual(graph.provides['/test_output.txt'].name, 'TestTarget')

    def test_raises_exceptions_if_two_targets_produce_the_same_file(self):
        self.workflow.target(
            'TestTarget1', inputs=[], outputs=['/test_output.txt'], working_dir='')
        self.workflow.target(
            'TestTarget2', inputs=[], outputs=['/test_output.txt'], working_dir='')

        with self.assertRaises(MultipleProvidersError):
            Graph.from_targets(self.workflow.targets)

    def test_finds_no_dependencies_for_target_with_no_inputs(self):
        target = self.workflow.target('TestTarget', inputs=[], outputs=[])
        graph = Graph.from_targets(self.workflow.targets)
        self.assertEqual(graph.dependencies[target], set())

    @patch('gwf.core.os.path.exists', return_value=False, autospec=True)
    def test_non_existing_files_not_provided_by_other_target_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget', inputs=['test_input.txt'], outputs=[])
        with self.assertRaises(MissingProviderError):
            Graph.from_targets(self.workflow.targets, )

    @patch('gwf.core.os.path.exists', return_value=True, autospec=True)
    def test_existing_files_not_provided_by_other_target_has_no_dependencies(self, mock_exists):
        target = self.workflow.target('TestTarget', inputs=['test_file.txt'], outputs=[])
        graph = Graph.from_targets(self.workflow.targets)
        self.assertEqual(graph.dependencies[target], set())

    @patch('gwf.core.os.path.exists', return_value=False, auto_spec=True)
    def test_finds_non_existing_file_provided_by_other_target(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=['test_file.txt'], outputs=[])

        graph = Graph.from_targets(self.workflow.targets)

        self.assertIn(target2, graph.dependencies)
        self.assertIn(target1, graph.dependencies[target2])

    @patch('gwf.core.os.path.exists', return_value=False, autospec=True)
    def test_finds_non_existing_files_provided_by_two_other_targets(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file1.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=[], outputs=['test_file2.txt'])
        target3 = self.workflow.target(
            'TestTarget3', inputs=['test_file1.txt', 'test_file2.txt'], outputs=[])

        graph = Graph.from_targets(self.workflow.targets)

        self.assertIn(target3, graph.dependencies)
        self.assertIn(target1, graph.dependencies[target3])
        self.assertIn(target2, graph.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False, autospec=True)
    def test_finds_non_existing_files_provided_by_other_targets_in_chain(self, mock_os_path_exists):
        target1 = self.workflow.target(
            'TestTarget1', inputs=[], outputs=['test_file1.txt'])
        target2 = self.workflow.target(
            'TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])
        target3 = self.workflow.target(
            'TestTarget3', inputs=['test_file2.txt'], outputs=[])

        graph = Graph.from_targets(self.workflow.targets)
        self.assertIn(target2, graph.dependencies)
        self.assertIn(target3, graph.dependencies)
        self.assertIn(target1, graph.dependencies[target2])
        self.assertIn(target2, graph.dependencies[target3])

    @patch('gwf.core.os.path.exists', return_value=False, autospec=True)
    def test_immediate_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target(
            'TestTarget1', inputs=['test_file2.txt'], outputs=['test_file1.txt'])
        self.workflow.target(
            'TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])

        with self.assertRaises(CircularDependencyError):
            Graph.from_targets(self.workflow.targets)

    @patch('gwf.core.os.path.exists', return_value=False, autospec=True)
    def test_circular_dependencies_in_workflow_raises_exception(self, mock_os_path_exists):
        self.workflow.target('TestTarget1', inputs=['test_file3.txt'], outputs=['test_file1.txt'])
        self.workflow.target('TestTarget2', inputs=['test_file1.txt'], outputs=['test_file2.txt'])
        self.workflow.target('TestTarget3', inputs=['test_file2.txt'], outputs=['test_file3.txt'])

        # TODO: Rewrite this test to not use from_targets to build the graph.
        with self.assertRaises(CircularDependencyError):
            Graph.from_targets(self.workflow.targets)

    def test_endpoints_only_includes_target_with_no_dependents(self):
        self.workflow.target('TestTarget1', inputs=[], outputs=['test.txt'])
        target2 = self.workflow.target('TestTarget2', inputs=['test.txt'], outputs=[])
        target3 = self.workflow.target('TestTarget3', inputs=['test.txt'], outputs=[])
        graph = Graph.from_targets(self.workflow.targets)
        self.assertSetEqual(graph.endpoints(), {target2, target3})

    def test_dependencies_correctly_resolved_for_named_workflow(self):
        workflow = Workflow(name='foo')
        target1 = workflow.target('TestTarget1', inputs=[], outputs=['test.txt'])
        target2 = workflow.target('TestTarget2', inputs=['test.txt'], outputs=[])

        other_workflow = Workflow(name='bar')
        other_workflow.include(workflow)
        other_target1 = other_workflow.target('TestTarget1', inputs=['test.txt'], outputs=[])

        graph = Graph.from_targets(other_workflow.targets)
        assert 'TestTarget1' in graph
        assert 'foo.TestTarget2' in graph
        assert 'foo.TestTarget2' in graph

    def test_raise_error_if_two_targets_in_different_namespaces_produce_the_same_file(self):
        w1 = Workflow(name='foo')
        w1.target('SayHello', inputs=[], outputs=['greeting.txt'])

        w2 = Workflow(name='bar')
        w2.target('SayHi', inputs=[], outputs=['greeting.txt'])
        w2.include(w1)

        with self.assertRaises(MultipleProvidersError):
            Graph.from_targets(w2.targets)


class TestShouldRun(unittest.TestCase):
    def setUp(self):
        workflow = Workflow(working_dir='/some/dir')
        self.target1 = workflow.target('TestTarget1', inputs=[], outputs=['test_output1.txt'])
        self.target2 = workflow.target('TestTarget2', inputs=['test_output1.txt'], outputs=['test_output2.txt'])
        self.target3 = workflow.target('TestTarget3', inputs=['test_output1.txt'], outputs=['test_output3.txt'])
        self.target4 = workflow.target('TestTarget4', inputs=['test_output2.txt', 'test_output3.txt'],
                                       outputs=['final_output.txt'])

        self.graph = Graph.from_targets(workflow.targets)
        self.backend = DummyBackend()
        self.scheduler = Scheduler(graph=self.graph, backend=self.backend)

    def test_target_should_run_if_one_of_its_dependencies_does_not_exist(self):
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(self.scheduler.should_run(self.target1))

        self.assertEqual(
            logs.output,
            [
                'DEBUG:gwf.core:TestTarget1 should run because one of its output files does not exist.'
            ]
        )

    def test_target_should_run_if_one_of_its_dependencies_should_run(self):
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(self.scheduler.should_run(self.target2))

        self.assertEqual(
            logs.output,
            [
                'DEBUG:gwf.core:TestTarget1 should run because one of its output files does not exist.',
                'DEBUG:gwf.core:TestTarget2 should run because one of its dependencies should run.'
            ]
        )

    def test_target_should_run_if_it_is_a_sink(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        graph = Graph.from_targets({'TestTarget': target})
        scheduler = Scheduler(graph=graph, backend=DummyBackend())
        with self.assertLogs(level='DEBUG') as logs:
            self.assertTrue(scheduler.schedule(target))
            self.assertEqual(logs.output, [
                'DEBUG:gwf.core:Scheduling target TestTarget.',
                'DEBUG:gwf.core:TestTarget should run because it is a sink.',
                'INFO:gwf.core:Submitting target TestTarget.',
            ])

    def test_target_should_not_run_if_it_is_a_source_and_all_outputs_exist(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target(
            'TestTarget1',
            inputs=[],
            outputs=['test_output1.txt', 'test_output2.txt']
        )

        graph = Graph.from_targets(workflow.targets)
        scheduler = Scheduler(graph=graph, backend=DummyBackend())

        mock_file_cache = {
            '/some/dir/test_output1.txt': 1,
            '/some/dir/test_output2.txt': 2
        }
        with patch.dict(scheduler._file_cache, mock_file_cache):
            self.assertFalse(scheduler.should_run(target))

    def test_should_run_if_any_input_file_is_newer_than_any_output_file(self):
        mock_file_cache = {
            '/some/dir/test_output1.txt': 0,
            '/some/dir/test_output2.txt': 1,
            '/some/dir/test_output3.txt': 3,
            '/some/dir/final_output.txt': 2,
        }

        with patch.dict(self.scheduler._file_cache, mock_file_cache):
            self.assertFalse(self.scheduler.should_run(self.target1))
            self.assertFalse(self.scheduler.should_run(self.target2))
            self.assertFalse(self.scheduler.should_run(self.target3))
            self.assertTrue(self.scheduler.should_run(self.target4))

    def test_should_run_not_run_if_all_outputs_are_newer_then_the_inputs(self):
        mock_file_cache = {
            '/some/dir/test_output1.txt': 0,
            '/some/dir/test_output2.txt': 1,
            '/some/dir/test_output3.txt': 3,
            '/some/dir/final_output.txt': 4,
        }

        with patch.dict(self.scheduler._file_cache, mock_file_cache):
            self.assertFalse(self.scheduler.should_run(self.target1))
            self.assertFalse(self.scheduler.should_run(self.target2))
            self.assertFalse(self.scheduler.should_run(self.target3))
            self.assertFalse(self.scheduler.should_run(self.target4))


@patch('gwf.core.os.path.exists', return_value=True, autospec=True)
def test_two_targets_producing_the_same_file_but_declared_with_rel_and_abs_path(mock_os_path_exists):
    workflow = Workflow(working_dir='/some/dir')
    workflow.target('TestTarget1', inputs=[], outputs=['/some/dir/test_output.txt'])
    workflow.target('TestTarget2', inputs=[], outputs=['test_output.txt'])

    with pytest.raises(MultipleProvidersError):
        Graph.from_targets(workflow.targets)


def test_scheduling_submitted_target(backend, monkeypatch):
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    graph = Graph.from_targets({'TestTarget': target})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)
    backend.submit(target, dependencies=set())
    assert len(backend.submit.call_args_list) == 1
    assert scheduler.schedule(target) == True
    assert len(backend.submit.call_args_list) == 1


def test_scheduling_unsubmitted_target(backend, monkeypatch):
    target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
    graph = Graph.from_targets({'TestTarget': target})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)
    assert scheduler.schedule(target) == True
    assert len(backend.submit.call_args_list) == 1
    assert call(target, dependencies=set()) in backend.submit.call_args_list


def test_scheduling_target_with_deps_that_are_not_submitted(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['test_output.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=['test_output.txt'], outputs=[], options={}, working_dir='/some/dir')
    graph = Graph.from_targets({'TestTarget1': target1, 'TestTarget2': target2})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)
    assert scheduler.schedule(target2) == True
    assert len(backend.submit.call_args_list) == 2
    assert call(target1, dependencies=set()) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list


def test_scheduling_target_with_deep_deps_that_are_not_submitted(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['test_output1.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=['test_output1.txt'], outputs=['test_output2.txt'], options={},
                     working_dir='/some/dir')
    target3 = Target('TestTarget3', inputs=['test_output2.txt'], outputs=['test_output3.txt'], options={},
                     working_dir='/some/dir')
    target4 = Target('TestTarget4', inputs=['test_output3.txt'], outputs=['final_output.txt'], options={},
                     working_dir='/some/dir')
    graph = Graph.from_targets({'target1': target1, 'target2': target2, 'target3': target3, 'target4': target4})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)
    assert scheduler.schedule(target4) == True
    assert len(backend.submit.call_args_list) == 4
    assert call(target1, dependencies=set()) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target2])) in backend.submit.call_args_list
    assert call(target4, dependencies=set([target3])) in backend.submit.call_args_list


def test_scheduling_branch_and_join_structure(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['output1.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=['output1.txt'], outputs=['output2.txt'], options={},
                     working_dir='/some/dir')
    target3 = Target('TestTarget3', inputs=['output1.txt'], outputs=['output3.txt'], options={},
                     working_dir='/some/dir')
    target4 = Target('TestTarget4', inputs=['output2.txt', 'output3.txt'], outputs=['final.txt'], options={},
                     working_dir='/some/dir')
    graph = Graph.from_targets({'target1': target1, 'target2': target2, 'target3': target3, 'target4': target4})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)
    assert scheduler.schedule(target4) == True
    assert len(backend.submit.call_args_list) == 4
    assert call(target1, dependencies=set([])) in backend.submit.call_args_list
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target4, dependencies=set([target3, target2])) in backend.submit.call_args_list


def test_scheduling_branch_and_join_structure_with_previously_submitted_dependency(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['output1.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=['output1.txt'], outputs=['output2.txt'], options={},
                     working_dir='/some/dir')
    target3 = Target('TestTarget3', inputs=['output1.txt'], outputs=['output3.txt'], options={},
                     working_dir='/some/dir')
    target4 = Target('TestTarget4', inputs=['output2.txt', 'output3.txt'], outputs=['final.txt'], options={},
                     working_dir='/some/dir')

    graph = Graph.from_targets({'target1': target1, 'target2': target2, 'target3': target3, 'target4': target4})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)

    backend.submit(target1, dependencies=set())

    assert scheduler.schedule(target4) == True
    assert len(backend.submit.call_args_list) == 4
    assert call(target2, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target4, dependencies=set([target3, target2])) in backend.submit.call_args_list


def test_scheduling_non_submitted_targets_that_should_not_run(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['test_output1.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=[], outputs=['test_output2.txt'], options={}, working_dir='/some/dir')
    target3 = Target('TestTarget3', inputs=['test_output1.txt', 'test_output2.txt'], outputs=['test_output3.txt'],
                     options={}, working_dir='/some/dir')
    graph = Graph.from_targets({'TestTarget1': target1, 'TestTarget2': target2, 'TestTarget3': target3})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: False)
    assert scheduler.schedule(target3) == False
    assert backend.submit.call_args_list == []


def test_scheduling_many_targets_calls_schedule_for_each_target(backend, monkeypatch):
    target1 = Target('TestTarget1', inputs=[], outputs=['test_output1.txt'], options={}, working_dir='/some/dir')
    target2 = Target('TestTarget2', inputs=[], outputs=['test_output2.txt'], options={}, working_dir='/some/dir')
    target3 = Target('TestTarget3', inputs=['test_output1.txt'], outputs=['test_output3.txt'], options={},
                     working_dir='/some/dir')
    target4 = Target('TestTarget4', inputs=['test_output2.txt'], outputs=['test_output4.txt'], options={},
                     working_dir='/some/dir')
    graph = Graph.from_targets(
        {'TestTarget1': target1, 'TestTarget2': target2, 'TestTarget3': target3, 'TestTarget4': target4})
    scheduler = Scheduler(graph=graph, backend=backend)
    monkeypatch.setattr(scheduler, 'should_run', lambda t: True)

    assert scheduler.schedule_many([target3, target4]) == [True, True]
    assert call(target4, dependencies=set([target2])) in backend.submit.call_args_list
    assert call(target3, dependencies=set([target1])) in backend.submit.call_args_list
    assert call(target2, dependencies=set()) in backend.submit.call_args_list
    assert call(target1, dependencies=set()) in backend.submit.call_args_list
