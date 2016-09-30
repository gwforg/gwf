from io import StringIO
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

        self.mock_print = self.create_patch(
            'builtins.print'
        )

        self.mock_backend = Mock(
            name='backend', spec_set=Backend
        )
        self.mock_backend.logs.return_value = StringIO('this is the log...')

        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints']
        )
        self.mock_workflow.targets = {
            'TestTarget': sentinel.TestTarget,
        }

        self.logs_command = LogsCommand()
        self.logs_command.backend = self.mock_backend
        self.logs_command.workflow = self.mock_workflow

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
        self.logs_command.config = {
            'target': 'TestTarget',
            'stderr': False,
        }

        self.logs_command.on_run()

        self.mock_backend.logs.assert_called_once_with(
            sentinel.TestTarget, stderr=False,
        )

        self.mock_print.assert_called_once_with(
            'this is the log...', end=''
        )

    def test_raises_exception_if_target_name_does_not_exist_in_workflow(self):
        self.logs_command.config = {
            'target': 'WrongTarget',
            'stderr': False,
        }

        with self.assertRaises(TargetDoesNotExistError) as ex:
            self.logs_command.on_run()
            self.assertEqual('WrongTarget', ex.target)

    def test_requests_stderr_log_when_stderr_option_is_true(self):
        self.logs_command.config = {
            'target': 'TestTarget',
            'stderr': True,
        }

        self.logs_command.on_run()

        self.mock_backend.logs.assert_called_once_with(
            sentinel.TestTarget, stderr=True,
        )

        self.mock_print.assert_called_once_with(
            'this is the log...', end=''
        )
