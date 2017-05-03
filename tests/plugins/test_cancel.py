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
        self.assertEqual(result.output, 'Cancelling target Target1.\nCancelling target Target2.\n')

    def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel'], input='y')
        self.assertEqual(result.output, 'Cancelling target Target1.\nCancelling target Target2.\n')

    def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(self):
        result = self.runner.invoke(main, ['-b', 'testing', 'cancel'], input='N')
        self.assertIn('Aborted!\n', result.output)
