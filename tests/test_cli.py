import logging
from unittest.mock import ANY, Mock, sentinel

from gwf.backends.testing import TestingBackend
from gwf.cli import App, main
from gwf.exceptions import GWFError
from gwf.ext import Extension, ExtensionManager

from . import GWFTestCase


class TestingPlugin(Extension):

    name = 'plugin1'

    def configure(self, *args, **kwargs):
        pass

    def setup_argument_parser(self, parser, subparsers):
        self.setup_subparser(
            subparsers,
            'foo',
            'Command for fooing.',
            self.handle_subcommand1
        )

        self.setup_subparser(
            subparsers,
            'boo',
            'Command for booing.',
            self.handle_subcommand2
        )

    def handle_subcommand1(self):
        pass

    def handle_subcommand2(self):
        pass


class TestApp(GWFTestCase):

    def setUp(self):
        self.testing_plugin1 = TestingPlugin()
        self.testing_plugin1.handle_subcommand1 = Mock()
        self.testing_plugin1.handle_subcommand2 = Mock()

        self.mock_plugins_manager = ExtensionManager('gwf.plugins')
        self.mock_plugins_manager.load_extensions = Mock()
        self.mock_plugins_manager.exts = {
            'plugin1': self.testing_plugin1,
        }

        self.mock_backends_manager = ExtensionManager('gwf.backends')
        self.mock_plugins_manager.load_extensions = Mock()
        self.mock_backends_manager.exts = {
            'backend1': TestingBackend(),
        }

        self.mock_sys_path = self.create_patch(
            'gwf.cli.sys.path'
        )

        self.mock_import_object = self.create_patch(
            'gwf.cli.import_object'
        )

        self.mock_atexit_register = self.create_patch(
            'gwf.cli.atexit.register'
        )

        self.app = App(
            plugins_manager=self.mock_plugins_manager,
            backends_manager=self.mock_backends_manager,
            config_files=[],
        )

        self.app.parser.print_help = Mock()
        self.app.parser.print_usage = Mock()

    def test_run_adds_workflow_directory_to_python_path_and_exits(self):
        with self.assertRaises(SystemExit) as e:
            self.app.run(['-f', '/some/dir/workflow.py', '-b', 'backend1'])
            self.mock_sys_path.insert.assert_called_once_with(0, '/some/dir')
            self.app.parser.print_help.assert_called_once_with()
            self.assertEqual(e.code, 1)

    def test_run_with_no_subcommand_prints_help_and_exits_with_code_1(self):
        with self.assertRaises(SystemExit) as e:
            self.app.run(['-f', '/some/dir/workflow.py', '-b', 'backend1'])
            self.app.parser.print_help.assert_called_once_with()
            self.assertEqual(e.code, 1)

    def test_run_with_subcommand_calls_subcommand_handler(self):
        self.app.run(['-b', 'backend1', 'foo'])
        self.testing_plugin1.handle_subcommand1.assert_called_once_with()

        self.app.run(['-b', 'backend1', 'boo'])
        self.testing_plugin1.handle_subcommand2.assert_called_once_with()

    def test_run_with_unknown_subcommand_prints_error(self):
        with self.assertRaises(SystemExit) as e:
            self.app.run(['-b', 'backend1', 'bar'])
            self.app.parser.print_usage.assert_called_once_with()
            self.assertEqual(e.code, 1)


class TestMain(GWFTestCase):

    def setUp(self):
        self.mock_colorama = self.create_patch(
            'gwf.cli.colorama'
        )

        self.mock_app = self.create_patch(
            'gwf.cli.App'
        )

        self.mock_extension_manager = self.create_patch(
            'gwf.cli.ExtensionManager'
        )

    def test_main_initializes_app(self):
        self.mock_extension_manager.side_effect = [
            sentinel.plugin_manager,
            sentinel.backend_manager,
        ]

        main()

        self.mock_colorama.init.assert_called_once_with()

        self.mock_app.assert_called_once_with(
            plugins_manager=sentinel.plugin_manager,
            backends_manager=sentinel.backend_manager,
            config_files=ANY,
        )

        self.mock_app.return_value.run.assert_called_once_with(ANY)

    def test_main_logs_any_raised_gwf_errors_and_exits(self):
        self.mock_app.side_effect = GWFError('this is the error')

        with self.assertLogs(level='ERROR') as log:
            with self.assertRaisesRegex(SystemExit, '1'):
                main()
            self.assertIn('this is the error', log.output[0])
