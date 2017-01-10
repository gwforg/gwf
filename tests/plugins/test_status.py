from unittest.mock import ANY, Mock, sentinel

from gwf.backends.base import Backend
from gwf.exceptions import TargetDoesNotExistError
from gwf.plugins.status import StatusCommand

from .. import GWFTestCase


class StatusCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.status.StatusCommand.setup_subparser'
        )

        self.status_command = StatusCommand()
        self.status_command.print_progress = lambda targets: None

        self.mock_backend = Mock(name='backend', spec_set=Backend)
        self.mock_get_active_backend = Mock(return_value=self.mock_backend)

        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints', 'dependencies'],
            dependencies=[]
        )
        self.mock_get_prepared_workflow = Mock(return_value=self.mock_workflow)

        self.mock_get_terminal_size = self.create_patch('gwf.plugins.status.os.get_terminal_size')
        self.mock_get_terminal_size.return_value = 80

    def test_sets_up_status_subcommand(self):
        self.status_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers, 'status', ANY, self.status_command.on_run,
        )

        mock_subparser = self.mock_setup_subparser.return_value
        mock_subparser.add_argument.assert_any_call(
            'targets', metavar="TARGET", nargs='*', help=ANY,
        )

    def test_on_run_runs_all_endpoints_if_no_targets_are_given(self):
        self.mock_workflow.endpoints.return_value = [
            sentinel.target1, sentinel.target2
        ]

        mock_config = {'targets': [], 'verbose': True}

        self.status_command.configure(
            get_active_backend=self.mock_get_active_backend,
            get_prepared_workflow=self.mock_get_prepared_workflow,
            config=mock_config,
        )

        self.status_command.on_run()

        self.mock_workflow.endpoints.assert_called_once_with()

    def test_on_run_runs_with_targets_given_in_config(self):
        self.mock_workflow.targets = {
            'target1': sentinel.target1,
            'target2': sentinel.target2,
        }

        mock_config = {
            'targets': ['target1', 'target2'],
        }

        self.status_command.configure(
            get_active_backend=self.mock_get_active_backend,
            get_prepared_workflow=self.mock_get_prepared_workflow,
            config=mock_config,
        )

        self.status_command.on_run()

        self.mock_workflow.endpoints.assert_not_called()

    def test_on_run_raises_exception_if_target_does_not_exist_in_workflow(self):
        self.mock_workflow.targets = {
            'target1': sentinel.target1,
        }

        mock_config = {
            'targets': ['target1', 'target2'],
        }

        self.status_command.configure(
            get_active_backend=self.mock_get_active_backend,
            get_prepared_workflow=self.mock_get_prepared_workflow,
            config=mock_config,
        )

        with self.assertRaises(TargetDoesNotExistError) as ex:
            self.status_command.on_run()
            self.assertEqual(ex.name, 'target2')
