import unittest
from unittest.mock import patch, create_autospec

from gwf.core import Workflow
from gwf.backends.base import Backend
from gwf.exceptions import GWFException
from gwf.scheduler import Scheduler


class TestScheduler(unittest.TestCase):

    def setUp(self):
        self.mock_backend = create_autospec(Backend)
        self.workflow = Workflow()

    def test_finds_no_providers_in_empty_workflow(self):
        scheduler = Scheduler(self.workflow, self.mock_backend)
        self.assertDictEqual(scheduler.provides, {})

    def test_finds_no_providers_in_workflow_with_no_producers(self):
        self.workflow.target('TestTarget', inputs=[], outputs=[])

        scheduler = Scheduler(self.workflow, self.mock_backend)
        self.assertDictEqual(scheduler.provides, {})

    def test_finds_provider_in_workflow_with_one_producer(self):
        self.workflow.target('TestTarget', inputs=[], outputs=['/test_output.txt'], working_dir='')

        scheduler = Scheduler(self.workflow, self.mock_backend)
        self.assertIn('/test_output.txt', scheduler.provides)
        self.assertEqual(scheduler.provides['/test_output.txt'].name, 'TestTarget')

    def test_finds_no_dependencies_for_target_with_no_inputs(self):
        target = self.workflow.target('TestTarget', inputs=[], outputs=[])
        scheduler = Scheduler(self.workflow, self.mock_backend)

        self.assertNotIn(target, scheduler.dependencies)

    @patch('gwf.scheduler.os.path.exists', return_value=True)
    def test_existing_files_cause_no_dependencies(self, mock_os_path_exists):
        target = self.workflow.target('TestTarget', inputs=['test_input.txt'], outputs=[])

        scheduler = Scheduler(self.workflow, self.mock_backend)
        self.assertNotIn(target, scheduler.dependencies)

    @patch('gwf.scheduler.os.path.exists', return_value=False)
    def test_non_existing_files_not_provided_by_other_target_raises_exception(self, mock_os_path_exists):
        self.workflow.target('TestTarget', inputs=['test_input.txt'], outputs=[])
        with self.assertRaises(GWFException):
            Scheduler(self.workflow, self.mock_backend)

    @patch('gwf.scheduler.os.path.exists', return_value=False)
    def test_find_non_existing_file_provided_by_other_target(self, mock_os_path_exists):
        target1 = self.workflow.target('TestTarget1', inputs=[], outputs=['test_file.txt'])
        target2 = self.workflow.target('TestTarget2', inputs=['test_file.txt'], outputs=[])
        
        scheduler = Scheduler(self.workflow, self.mock_backend)

        self.assertIn(target2, scheduler.dependencies)
        self.assertIn(target1, scheduler.dependencies[target2])

    @patch('gwf.scheduler.os.path.exists', return_value=False)
    def test_find_non_existing_files_provided_by_two_other_targets(self, mock_os_path_exists):
        target1 = self.workflow.target('TestTarget1', inputs=[], outputs=['test_file1.txt'])
        target2 = self.workflow.target('TestTarget2', inputs=[], outputs=['test_file2.txt'])
        target3 = self.workflow.target('TestTarget3', inputs=['test_file1.txt', 'test_file2.txt'], outputs=[])

        scheduler = Scheduler(self.workflow, self.mock_backend)

        self.assertIn(target3, scheduler.dependencies)
        self.assertIn(target1, scheduler.dependencies[target3])
        self.assertIn(target2, scheduler.dependencies[target3])

    def test_raises_exceptions_if_two_targets_produce_the_same_file(self):
        self.workflow.target('TestTarget1', inputs=[], outputs=['/test_output.txt'], working_dir='')
        self.workflow.target('TestTarget2', inputs=[], outputs=['/test_output.txt'], working_dir='')

        with self.assertRaises(GWFException):
            Scheduler(self.workflow, self.mock_backend)
