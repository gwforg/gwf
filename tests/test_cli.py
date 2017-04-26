from tests import CliTestCase
from gwf.cli import main


class TestCli(CliTestCase):

    def test_main_shows_usage_when_no_subcommand_is_given(self):
        result = self.runner.invoke(main, [])
        self.assertEqual(result.exit_code, 0)
        self.assertTrue(result.output.startswith('Usage:'))