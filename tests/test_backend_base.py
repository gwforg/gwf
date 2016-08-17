import unittest

from gwf.backends.base import Backend, get_backends
from gwf.core import Workflow
from gwf.exceptions import WorkflowNotPreparedError


class TestBaseBackend(unittest.TestCase):

    def test_subclass_of_Backend_is_registered_in_backend_registry(self):
        class TestingBackend(Backend):
            pass

        backends = get_backends()
        self.assertIn('TestingBackend', backends)
        self.assertEqual(TestingBackend, backends['TestingBackend'])

    def test_backend_raises_exception_when_initialized_with_unprepared_workflow(self):
        with self.assertRaises(WorkflowNotPreparedError):
            Backend(workflow=Workflow())
