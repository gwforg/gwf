from unittest.mock import ANY, call, mock_open, patch, sentinel

from gwf.exceptions import GWFError
from gwf.plugins.config import ConfigCommand

from .. import GWFTestCase


class ConfigCommandTest(GWFTestCase):

    def setUp(self):
        self.mock_setup_subparser = self.create_patch(
            'gwf.plugins.config.ConfigCommand.setup_subparser'
        )
        self.config_command = ConfigCommand()

    def test_setup_argument_parser_correctly_sets_up_subcommand(self):
        self.config_command.setup_argument_parser(
            sentinel.parser, sentinel.subparsers,
        )

        self.mock_setup_subparser.assert_called_once_with(
            sentinel.subparsers,
            'config',
            ANY,
            self.config_command.on_run,
        )

        parser = self.mock_setup_subparser.return_value
        parser.add_argument.assert_has_calls([
            call('-u', '--user', action='store_true', help=ANY),
            call('option_name', metavar='option', help=ANY),
            call('option_value', metavar='value',
                 nargs='?', default='', help=ANY),
        ])

    def test_set_an_option_on_non_existing_local_conf_file(self):
        self.config_command.config = {
            'option_name': 'foo',
            'option_value': 'bar',
            'user': False
        }

        m = mock_open()
        with patch('builtins.open', m) as mock_open_:
            mock_open_.side_effect = [IOError, m.return_value]

            self.config_command.on_run()

            self.assertIn(call('.gwf.conf'), mock_open_.call_args_list)
            self.assertIn(call('.gwf.conf', 'w'), mock_open_.call_args_list)

            self.assertIn(
                call('foo = bar\n'),
                mock_open_.return_value.write.call_args_list
            )

    def test_set_an_option_on_existing_local_conf_file(self):
        self.config_command.config = {
            'option_name': 'foo',
            'option_value': 'bar',
            'user': False
        }

        with patch('builtins.open', mock_open(read_data='boz = jax\n')) as mock_open_:
            self.config_command.on_run()

            self.assertIn(call('.gwf.conf'), mock_open_.call_args_list)
            self.assertIn(call('.gwf.conf', 'w'), mock_open_.call_args_list)

            self.assertIn(
                call('boz = jax\nfoo = bar\n'),
                mock_open_.return_value.write.call_args_list
            )

    def test_unset_existing_option(self):
        self.config_command.config = {
            'option_name': 'foo',
            'option_value': '',
            'user': False
        }

        with patch('builtins.open', mock_open(read_data='foo = bar\n')) as mock_open_:
            self.config_command.on_run()

            self.assertIn(
                call(''),
                mock_open_.return_value.write.call_args_list
            )

    def test_unset_non_existing_option(self):
        self.config_command.config = {
            'option_name': 'koo',
            'option_value': '',
            'user': False
        }

        with patch('builtins.open', mock_open(read_data='foo = bar\n')):
            with self.assertRaises(GWFError):
                self.config_command.on_run()
