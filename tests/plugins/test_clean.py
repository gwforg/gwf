from unittest.mock import ANY, Mock, call, patch, sentinel

from gwf.backends.base import Backend
from gwf.exceptions import TargetDoesNotExistError
from gwf.plugins.clean import CleanCommand, _delete_file

from .. import GWFTestCase


class CleanCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.clean.CleanCommand.setup_subparser'
        )

        self.mock_backend = Mock(name='backend', spec_set=Backend)
        self.mock_get_active_backend = Mock(return_value=self.mock_backend)

        self.clean_command = CleanCommand()

        self.mock_workflow = Mock(
            name='workflow',
            spec_set=['targets', 'endpoints', 'working_dir']
        )
        self.mock_get_prepared_workflow = Mock(return_value=self.mock_workflow)

        self.mock_target1 = Mock(spec_set=['name', 'outputs'])
        self.mock_target1.name = 'target1'
        self.mock_target1.outputs = ['/some/dir/foo.txt']

        self.mock_target2 = Mock(spec_set=['name', 'outputs'])
        self.mock_target2.name = 'target2'
        self.mock_target2.outputs = ['/some/dir/bar.txt']

        self.mock_workflow.targets = {
            'target1': self.mock_target1,
            'target2': self.mock_target2,
        }
        self.mock_workflow.working_dir = '/some/dir'

        self.mock_delete_file = self.create_patch(
            'gwf.plugins.clean._delete_file'
        )

    def test_sets_up_clean_subcommand(self):
        self.clean_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers, 'clean', ANY, self.clean_command.on_clean,
        )

        mock_subparser = self.mock_setup_subparser.return_value
        mock_subparser.add_argument.assert_has_calls([
            call('targets', metavar=ANY, nargs='*', help=ANY),
            call('-f', '--only-failed', action=ANY, help=ANY)
        ])

    def test_on_clean_cleans_all_targets_if_no_targets_are_given(self):
        mock_config = {'targets': [], 'only_failed': False}

        self.clean_command.configure(
            get_prepared_workflow=self.mock_get_prepared_workflow,
            get_active_backend=self.mock_get_active_backend,
            config=mock_config
        )

        self.clean_command.on_clean()

        self.assertIn(call('/some/dir/bar.txt'),
                      self.mock_delete_file.call_args_list)
        self.assertIn(call('/some/dir/foo.txt'),
                      self.mock_delete_file.call_args_list)

    def test_on_clean_cleans_with_targets_given_in_config(self):
        mock_config = {'targets': ['target1'], 'only_failed': False}

        self.clean_command.configure(
            get_prepared_workflow=self.mock_get_prepared_workflow,
            get_active_backend=self.mock_get_active_backend,
            config=mock_config
        )

        self.clean_command.on_clean()

        self.mock_delete_file.assert_has_calls([
            call('/some/dir/foo.txt'),
        ])

    def test_on_clean_raises_exception_if_target_does_not_exist_in_workflow(self):
        mock_config = {'targets': ['target1', 'target3'], 'only_failed': False}

        self.clean_command.configure(
            get_prepared_workflow=self.mock_get_prepared_workflow,
            get_active_backend=self.mock_get_active_backend,
            config=mock_config
        )

        with self.assertRaises(TargetDoesNotExistError) as ex:
            self.clean_command.on_clean()
            self.assertEqual(ex.name, 'target3')

    def test_on_clean_with_only_failed_flag_only_cleans_failed_targets(self):
        mock_config = {'targets': ['target1', 'target2'], 'only_failed': True}
        self.mock_backend.failed.side_effect = [False, True]

        self.clean_command.configure(
            get_prepared_workflow=self.mock_get_prepared_workflow,
            get_active_backend=self.mock_get_active_backend,
            config=mock_config
        )
        self.clean_command.on_clean()

        self.assertIn(call('/some/dir/bar.txt'),
                      self.mock_delete_file.call_args_list)
        self.assertNotIn(call('/some/dir/foo.txt'),
                         self.mock_delete_file.call_args_list)

    @patch('gwf.plugins.clean.os.remove', side_effect=IOError)
    def test_delete_file_ignores_non_existing_file(self, mock_os_remove):
        _delete_file('/some/dir/foo.txt')
        mock_os_remove.assert_called_once_with('/some/dir/foo.txt')

    @patch('gwf.plugins.clean.os.remove')
    def test_delete_file_deletes_existing_file(self, mock_os_remove):
        _delete_file('/some/dir/foo.txt')
        mock_os_remove.assert_called_once_with('/some/dir/foo.txt')
