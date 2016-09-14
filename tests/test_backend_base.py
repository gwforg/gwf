import logging
import unittest
from unittest.mock import patch

from gwf.backends.testing import TestingBackend
from gwf.core import PreparedWorkflow, Workflow
from gwf.exceptions import WorkflowNotPreparedError


class TestBaseBackend(unittest.TestCase):

    def test_backend_raises_exception_when_initialized_with_unprepared_workflow(self):
        with self.assertRaises(WorkflowNotPreparedError):
            TestingBackend(workflow=Workflow())

    def test_scheduling_workflow_with_one_target_that_is_already_submitted_returns_empty_schedule(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        backend = TestingBackend(PreparedWorkflow(workflow))
        with patch.object(backend, 'submitted', return_value=True):
            schedule = backend.schedule(target)
            self.assertListEqual(schedule, [])

    def test_scheduling_workflow_with_one_target_that_is_not_submitted_returns_schedule_with_target(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target)
                self.assertListEqual(schedule, [target])

    def test_scheduling_workflow_with_target_with_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1', outputs=['test_output.txt'])
        target2 = workflow.target('TestTarget2', inputs=['test_output.txt'])

        prepared_workflow = PreparedWorkflow(workflow)
        backend = TestingBackend(prepared_workflow)
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
        backend = TestingBackend(prepared_workflow)
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
        backend = TestingBackend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule = backend.schedule(target4)
                self.assertListEqual(
                    schedule, [target1, target2, target3, target4])

    def test_warn_if_target_uses_option_not_supported_by_backend(self):
        workflow = Workflow(working_dir='/some/dir')
        workflow.target('TestTarget', inputs=[],
                        outputs=[], walltime='01:00:00')

        prepared_workflow = PreparedWorkflow(workflow)
        with self.assertLogs(level=logging.WARN) as logs:
            TestingBackend(prepared_workflow)

            self.assertIn(
                'Backend "testing" does not support option "walltime".',
                logs.output[0]
            )
