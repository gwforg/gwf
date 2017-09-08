from tests import CliTestCase, touch_file

from gwf.cli import main

SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << "echo hello world"
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt']) << "echo world hello"
"""


class TestStatus(CliTestCase):

    def setUp(self):
        super().setUp()
        touch_file('workflow.py', contents=SIMPLE_WORKFLOW)

    def test_status_only_shows_endpoints_by_default(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'status'])
        self.assertIn('Target2', result.output)
        self.assertNotIn('Target1', result.output)

    def test_status_shows_all_targets_when_all_flag_is_used(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'status', '--all'])
        self.assertIn('Target2', result.output)
        self.assertIn('Target1', result.output)

    def test_status_shows_one_named_target(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'status', 'Target1'])
        self.assertNotIn('Target2', result.output)
        self.assertIn('Target1', result.output)

    def test_status_shows_two_named_targets(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'status', 'Target1', 'Target2'])
        self.assertIn('Target2', result.output)
        self.assertIn('Target1', result.output)

    def test_status_only_shows_names_when_only_names_flag_is_used(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'status', '--names-only'])
        self.assertEqual(result.output, 'Target2\n')