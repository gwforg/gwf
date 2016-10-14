from unittest.mock import ANY, Mock, sentinel

from gwf.backends.base import Backend
from gwf.exceptions import TargetDoesNotExistError
from gwf.plugins.run import RunCommand

from .. import GWFTestCase


class RunCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.run.RunCommand.setup_subparser'
        )

        self.mock_schedule_many = self.create_patch(
            'gwf.plugins.run.schedule_many'
        )

        self.run_command = RunCommand()

        self.mock_backend = Mock(name='backend', spec_set=Backend)
        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints']
        )

    def test_sets_up_run_subcommand(self):
        self.run_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers, 'run', ANY, self.run_command.on_run,
        )

        mock_subparser = self.mock_setup_subparser.return_value
        mock_subparser.add_argument.assert_called_once_with(
            'targets', metavar=ANY, nargs='*', help=ANY,
        )

    def test_on_run_runs_all_endpoints_if_no_targets_are_given(self):
        self.mock_workflow.endpoints.return_value = [
            sentinel.target1, sentinel.target2
        ]

        mock_config = {'targets': []}

        self.run_command.configure(workflow=self.mock_workflow,
                                   backend=self.mock_backend,
                                   config=mock_config)

        self.run_command.on_run()

        self.mock_workflow.endpoints.assert_called_once_with()
        self.mock_schedule_many.assert_called_once_with(
            self.mock_workflow,
            self.mock_backend,
            [sentinel.target1, sentinel.target2]
        )

    def test_on_run_runs_with_targets_given_in_config(self):
        self.mock_workflow.targets = {
            'target1': sentinel.target1,
            'target2': sentinel.target2,
        }

        mock_config = {'targets': ['target1', 'target2']}

        self.run_command.configure(workflow=self.mock_workflow,
                                   backend=self.mock_backend,
                                   config=mock_config)

        self.run_command.on_run()

        self.mock_workflow.endpoints.assert_not_called()
        self.mock_schedule_many.assert_called_once_with(
            self.mock_workflow, self.mock_backend,
            [sentinel.target1, sentinel.target2]
        )

    def test_on_run_raises_exception_if_target_does_not_exist_in_workflow(self):
        self.mock_workflow.targets = {
            'target1': sentinel.target1,
        }

        mock_config = {'targets': ['target1', 'target2']}

        self.run_command.configure(workflow=self.mock_workflow,
                                   backend=self.mock_backend,
                                   config=mock_config)

        with self.assertRaises(TargetDoesNotExistError) as ex:
            self.run_command.on_run()
            self.assertEqual(ex.name, 'target2')
