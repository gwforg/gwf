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
        self.mock_plugins_manager.close_extensions = Mock()
        self.mock_plugins_manager.exts = {
            'plugin1': self.testing_plugin1,
        }

        self.testing_backend1 = TestingBackend(working_dir='/some/dir')
        self.testing_backend1.close = Mock()

        self.mock_backends_manager = ExtensionManager('gwf.backends')
        self.mock_plugins_manager.load_extensions = Mock()
        self.mock_backends_manager.exts = {
            'backend1': self.testing_backend1,
        }

        self.mock_import_object = self.create_patch(
            'gwf.cli.import_object'
        )

        self.mock_ensure_dir = self.create_patch(
            'gwf.cli.ensure_dir'
        )

        self.app = App(
            plugins_manager=self.mock_plugins_manager,
            backends_manager=self.mock_backends_manager,
            config_files=[],
        )

        self.app.parser.print_help = Mock()
        self.app.parser.print_usage = Mock()


class TestMain(GWFTestCase):

    def setUp(self):
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

        self.mock_app.assert_called_once_with(
            plugins_manager=sentinel.plugin_manager,
            backends_manager=sentinel.backend_manager,
        )

        self.mock_app.return_value.run.assert_called_once_with(ANY)

    def test_main_logs_any_raised_gwf_errors_and_exits(self):
        self.mock_app.return_value.run.side_effect = GWFError('an error')

        with self.assertLogs(level='ERROR') as log:
            main()

        self.assertIn('an error', log.output[0])
        self.mock_app.return_value.exit.assert_called_once_with()
