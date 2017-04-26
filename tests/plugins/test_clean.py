import os.path

from tests import CliTestCase, touch_file
from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << "echo hello world"
gwf.target('Target2', inputs=[], outputs=['b.txt']) << "echo world hello"
"""


class TestRun(CliTestCase):

    def setUp(self):
        super().setUp()

        touch_file('workflow.py', contents=SIMPLE_WORKFLOW)
        touch_file('a.txt')
        touch_file('b.txt')

    def test_clean_output_from_all_targets_by_default(self):
        args = ['-b', 'testing', 'clean']
        self.runner.invoke(main, args)

        self.assertFalse(os.path.exists('a.txt'))
        self.assertFalse(os.path.exists('b.txt'))

    def test_clean_output_from_single_target(self):
        args = ['-b', 'testing', 'clean', 'Target1']
        self.runner.invoke(main, args)

        self.assertFalse(os.path.exists('a.txt'))
        self.assertTrue(os.path.exists('b.txt'))

    def test_clean_output_from_two_targets(self):
        args = ['-b', 'testing', 'clean', 'Target1', 'Target2']
        self.runner.invoke(main, args)

        self.assertFalse(os.path.exists('a.txt'))
        self.assertFalse(os.path.exists('b.txt'))

    def test_do_not_clean_outputs_from_endpoints(self):
        args = ['-b', 'testing', 'clean', '--not-endpoints', 'Target1', 'Target2']
        self.runner.invoke(main, args)

        self.assertTrue(os.path.exists('a.txt'))
        self.assertTrue(os.path.exists('b.txt'))
