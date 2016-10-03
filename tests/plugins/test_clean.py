from unittest.mock import ANY, Mock, sentinel

from gwf.backends.base import Backend
from gwf.exceptions import TargetDoesNotExistError
from gwf.plugins.clean import CleanCommand

from .. import GWFTestCase


class CleanCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.clean.CleanCommand.setup_subparser'
        )

        self.mock_backend = Mock(name='backend', spec_set=Backend)

        self.clean_command = CleanCommand()

        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints']
        )

        self.mock_target1 = Mock(spec_set=['clean', 'name'])
        self.mock_target1.name = 'target1'

        self.mock_target2 = Mock(spec_set=['clean', 'name'])
        self.mock_target2.name = 'target2'

        self.mock_workflow.targets = {
            'target1': self.mock_target1,
            'target2': self.mock_target2,
        }

    def test_sets_up_clean_subcommand(self):
        self.clean_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers, 'clean', ANY, self.clean_command.on_clean,
        )

        mock_subparser = self.mock_setup_subparser.return_value
        mock_subparser.add_argument.assert_called_once_with(
            'targets', metavar=ANY, nargs='*', help=ANY,
        )

    def test_on_clean_cleans_all_targets_if_no_targets_are_given(self):
        mock_config = {'targets': []}

        self.clean_command.configure(workflow=self.mock_workflow,
                                     backend=self.mock_backend,
                                     config=mock_config)

        self.clean_command.on_clean()

        self.mock_target1.clean.assert_called_once_with()
        self.mock_target2.clean.assert_called_once_with()

    def test_on_clean_cleans_with_targets_given_in_config(self):
        mock_config = {'targets': ['target1']}

        self.clean_command.configure(workflow=self.mock_workflow,
                                     backend=self.mock_backend,
                                     config=mock_config)

        self.clean_command.on_clean()

        self.mock_target1.clean.assert_called_once_with()
        self.mock_target2.clean.assert_not_called()

    def test_on_clean_raises_exception_if_target_does_not_exist_in_workflow(self):
        mock_config = {'targets': ['target1', 'target3']}

        self.clean_command.configure(workflow=self.mock_workflow,
                                     backend=self.mock_backend,
                                     config=mock_config)

        with self.assertRaises(TargetDoesNotExistError) as ex:
            self.clean_command.on_clean()
            self.assertEqual(ex.name, 'target3')
