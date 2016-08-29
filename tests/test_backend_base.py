import unittest
from unittest.mock import patch

from gwf.backends.base import Backend, get_backends
from gwf.core import PreparedWorkflow, Workflow
from gwf.exceptions import GWFError, WorkflowNotPreparedError


class TestBaseBackend(unittest.TestCase):

    def test_raises_exception_when_backend_does_not_declare_name(self):
        with self.assertRaises(GWFError):
            class TestingBackend(Backend):
                pass

    def test_subclass_of_Backend_is_registered_in_backend_registry(self):
        class TestingBackend(Backend):
            name = 'testing'

        backends = get_backends()
        self.assertIn('testing', backends)
        self.assertEqual(TestingBackend, backends['testing'])

    def test_backend_raises_exception_when_initialized_with_unprepared_workflow(self):
        with self.assertRaises(WorkflowNotPreparedError):
            Backend(workflow=Workflow())

    def test_scheduling_workflow_with_one_target_that_is_already_submitted_returns_empty_schedule(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        backend = Backend(PreparedWorkflow(workflow))
        with patch.object(backend, 'submitted', return_value=True):
            schedule, scheduled = backend.schedule(target)
            self.assertListEqual(schedule, [])
            self.assertSetEqual(scheduled, set())

    def test_scheduling_workflow_with_one_target_that_is_not_submitted_returns_schedule_with_target(self):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)
        backend = Backend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule, scheduled = backend.schedule(target)
                self.assertListEqual(schedule, [target])
                self.assertSetEqual(scheduled, set([target]))

    def test_scheduling_workflow_with_target_with_deps_that_are_not_submitted(self):
        workflow = Workflow(working_dir='/some/dir')
        target1 = workflow.target('TestTarget1', outputs=['test_output.txt'])
        target2 = workflow.target('TestTarget2', inputs=['test_output.txt'])

        prepared_workflow = PreparedWorkflow(workflow)
        backend = Backend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule, scheduled = backend.schedule(target2)
                self.assertListEqual(schedule, [target1, target2])
                self.assertSetEqual(scheduled, set([target1, target2]))

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
        backend = Backend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule, scheduled = backend.schedule(target4)
                self.assertListEqual(
                    schedule, [target1, target2, target3, target4])
                self.assertSetEqual(scheduled, set(
                    [target1, target2, target3, target4]))

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
        backend = Backend(prepared_workflow)
        with patch.object(backend, 'submitted', return_value=False):
            with patch.object(prepared_workflow, 'should_run', return_value=True):
                schedule, scheduled = backend.schedule(target4)
                self.assertListEqual(
                    schedule, [target1, target2, target3, target4])
                self.assertSetEqual(scheduled, set(
                    [target1, target2, target3, target4]))
