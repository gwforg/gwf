import logging
import unittest
from unittest.mock import patch

from gwf.backends.testing import TestingBackend
from gwf.core import PreparedWorkflow, Workflow


class TestBaseBackend(unittest.TestCase):

    def test_scheduling_workflow_with_one_target_that_is_already_submitted_returns_empty_schedule(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        backend = TestingBackend()
        backend.configure(workflow=PreparedWorkflow(workflow), config={})
        with patch.object(backend, 'submitted', return_value=True):
            schedule = backend.schedule(target)
            self.assertListEqual(schedule, [])

    def test_scheduling_workflow_with_one_target_that_is_not_submitted_returns_schedule_with_target(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target)
                self.assertListEqual(schedule, [target])

    def test_scheduling_workflow_with_target_with_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1', outputs=['test_output.txt'])
        target2 = workflow.target('TestTarget2', inputs=['test_output.txt'])

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target2)
                self.assertListEqual(schedule, [target1, target2])

    def test_scheduling_workflow_with_deep_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1', outputs=['test_output1.txt'])
        target2 = workflow.target('TestTarget2', inputs=[
                                  'test_output1.txt'], outputs=['test_output2.txt'])
        target3 = workflow.target('TestTarget3', inputs=[
                                  'test_output2.txt'], outputs=['test_output3.txt'])
        target4 = workflow.target('TestTarget4', inputs=[
                                  'test_output3.txt'], outputs=['final_output.txt'])

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target4)
                self.assertListEqual(
                    schedule, [target1, target2, target3, target4])

    def test_scheduling_workflow_with_branch_and_join_structure(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1', outputs=['test_output1.txt'])
        target2 = workflow.target('TestTarget2', inputs=[
                                  'test_output1.txt'], outputs=['test_output2.txt'])
        target3 = workflow.target('TestTarget3', inputs=[
                                  'test_output1.txt'], outputs=['test_output3.txt'])
        target4 = workflow.target(
            'TestTarget4',
            inputs=['test_output2.txt', 'test_output3.txt'],
            outputs=['final_output.txt']
        )

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target4)
                self.assertListEqual(
                    schedule, [target1, target2, target3, target4])

    @unittest.skip('Reimplement this check again.')
    def test_warn_if_target_uses_option_not_supported_by_backend(self):
        workflow = Workflow(working_dir='/some/dir')
        workflow.target('TestTarget', inputs=[],
                        outputs=[], walltime='01:00:00')

        prepared_workflow = PreparedWorkflow(workflow)
        with self.assertLogs(level=logging.WARN) as logs:
            backend = TestingBackend()
            backend.configure(workflow=prepared_workflow, config={})

            self.assertIn(
                'Backend "testing" does not support option "walltime".',
                logs.output[0]
            )

    def test_option_defaults_contains_specified_option_defaults(self):
        workflow = Workflow(working_dir='/some/dir')
        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        self.assertDictEqual(
            backend.option_defaults,
            {'cores': 2, 'memory': '18gb'}
        )

    def test_scheduling_a_submitted_dependency_does_not_submit_the_dependency(self):
        workflow = Workflow(working_dir='/some/dir')

        target1 = workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt']
        )
        target2 = workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        )
        target3 = workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        )

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', side_effect=[False, True, False, False]):
            schedule = backend.schedule(target3)
            self.assertEqual(schedule, [target2, target3])

    def test_scheduling_non_submitted_targets_that_should_not_run_does_not_submit_any_targets(self):
        workflow = Workflow(working_dir='/some/dir')

        target1 = workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt']
        )
        target2 = workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        )
        target3 = workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        )

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', side_effect=[False, False, False, False]):
                schedule = backend.schedule(target3)
                self.assertEqual(schedule, [])

    def test_scheduling_many_targets_calls_schedule_for_each_target(self):
        workflow = Workflow(working_dir='/some/dir')

        target1 = workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt']
        )
        target2 = workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        )
        target3 = workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt'],
            outputs=['test_output3.txt']
        )
        target4 = workflow.target(
            'TestTarget4',
            inputs=['test_output2.txt'],
            outputs=['test_output4.txt']
        )

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend()
        backend.configure(workflow=prepared_workflow, config={})
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule_many([target3, target4])
                self.assertEqual(
                    schedule, [[target1, target3], [target2, target4]])
