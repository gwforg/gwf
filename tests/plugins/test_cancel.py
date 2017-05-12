from tests import CliTestCase, touch_file
from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << "echo hello world"
gwf.target('Target2', inputs=[], outputs=['b.txt']) << "echo world hello"
"""


class TestCancel(CliTestCase):

    def setUp(self):
        super().setUp()
        touch_file('workflow.py', contents=SIMPLE_WORKFLOW)

    def test_cancel_one_target(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel', 'Target1'])
        self.assertEqual(result.output, 'Cancelling target Target1.\n')

    def test_cancel_two_targets(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel', 'Target1', 'Target2'])
        lines = result.output.split('\n')
        self.assertEqual(len(lines), 3)
        self.assertIn('Cancelling target Target1.', lines)
        self.assertIn('Cancelling target Target2.', lines)

    def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel'], input='y')
        lines = result.output.split('\n')
        self.assertEqual(len(lines), 3)
        self.assertIn('Cancelling target Target1.', lines)
        self.assertIn('Cancelling target Target2.', lines)

    def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel'], input='N')
        self.assertIn('Aborted!\n', result.output)
