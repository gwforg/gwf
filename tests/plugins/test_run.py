from unittest.mock import patch

from tests import CliTestCase, touch_file
from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=[]) << "echo hello world"
gwf.target('Target2', inputs=[], outputs=[]) << "echo world hello"
"""


@patch('gwf.plugins.run.schedule_many')
class TestRun(CliTestCase):

    def setUp(self):
        super().setUp()
        touch_file('workflow.py', SIMPLE_WORKFLOW)

    def test_run_all_targets(self, mock_schedule_many):
        args = ['-b', 'testing', 'run']
        self.runner.invoke(main, args)

        (graph, backend, targets), kwargs = mock_schedule_many.call_args
        self.assertEqual(len(targets), 2)
        self.assertEqual({x.name for x in targets}, {'Target1', 'Target2'})

    def test_run_specified_target(self, mock_schedule_many):
        args = ['-b', 'testing', 'run', 'Target1']
        self.runner.invoke(main, args)

        (graph, backend, targets), kwargs = mock_schedule_many.call_args
        self.assertEqual(len(targets), 1)
        self.assertEqual({x.name for x in targets}, {'Target1'})

    def test_dry_run(self, mock_schedule_many):
        args = ['-b', 'testing', 'run', '--dry-run']
        self.runner.invoke(main, args)

        (graph, backend, targets), kwargs = mock_schedule_many.call_args
        self.assertTrue(kwargs['dry_run'])