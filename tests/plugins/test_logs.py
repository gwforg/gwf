from unittest.mock import ANY, Mock, call, sentinel

from gwf.backends.base import Backend
from gwf.exceptions import TargetDoesNotExistError
from gwf.plugins.logs import LogsCommand

from .. import GWFTestCase


class LogsCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.logs.LogsCommand.setup_subparser'
        )

        self.logs_command = LogsCommand()

        self.mock_backend = Mock(name='backend', spec_set=Backend)
        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints']
        )

    def test_sets_up_logs_subcommand(self):
        self.logs_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers, 'logs', ANY, self.logs_command.on_run,
        )

        mock_subparser = self.mock_setup_subparser.return_value

        mock_subparser.add_argument.assert_has_calls([
            call('-e', '--stderr', action='store_true', help=ANY),
            call('target', help=ANY),
        ])

    def test_shows_stdout_log_for_target_in_workflow(self):
        pass
