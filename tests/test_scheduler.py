import unittest
from unittest.mock import patch, create_autospec

from gwf.core import Workflow
from gwf.backends.base import Backend
from gwf.scheduler import Scheduler


class TestScheduler(unittest.TestCase):

    def setUp(self):
        self.mock_backend = create_autospec(Backend)

    def test_scheduler_finds_no_providers_in_empty_workflow(self):
        workflow = Workflow()
        scheduler = Scheduler(workflow, self.mock_backend)
        self.assertDictEqual(scheduler.provides, {})

    def test_scheduler_finds_no_providers_in_workflow_with_no_producers(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=[])

        scheduler = Scheduler(workflow, self.mock_backend)
        self.assertDictEqual(scheduler.provides, {})

    def test_scheduler_finds_provider_in_workflow_with_one_producer(self):
        workflow = Workflow()
        workflow.target('TestTarget', inputs=[], outputs=['/test_output.txt'], working_dir='')

        scheduler = Scheduler(workflow, self.mock_backend)
        self.assertIn('/test_output.txt', scheduler.provides)
        self.assertEqual(scheduler.provides['/test_output.txt'].name, 'TestTarget')
