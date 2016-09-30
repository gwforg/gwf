from unittest.mock import Mock, sentinel

from gwf.exceptions import GWFError
from gwf.ext import ExtensionManager

from . import GWFTestCase


class TestExtensionManager(GWFTestCase):

    def setUp(self):
        self.mock_iter_entry_points = self.create_patch(
            'gwf.ext.iter_entry_points'
        )

        self.plugin1 = Mock(name='plugin1', spec_set=['name', '__init__'])
        self.plugin1.name = 'plugin1'
        self.plugin1.return_value = Mock(name='plugin1_inst')

        self.plugin2 = Mock(name='plugin2', spec_set=['name', '__init__'])
        self.plugin2.name = 'plugin2'
        self.plugin1.return_value = Mock(name='plugin2_inst')

        self.entrypoint1 = Mock(name='entrypoint1', spec_set=['load'])
        self.entrypoint1.load.return_value = self.plugin1

        self.entrypoint2 = Mock(name='entrypoint2', spec_set=['load'])
        self.entrypoint2.load.return_value = self.plugin2

        self.mock_iter_entry_points.return_value = [
            self.entrypoint1,
            self.entrypoint2,
        ]

    def test_load_extensions_loads_and_initializes_all_entry_points_in_group(self):
        em = ExtensionManager(group='myplugins')
        em.load_extensions()

        self.mock_iter_entry_points.assert_called_once_with(
            group='myplugins',
            name=None,
        )

        self.entrypoint1.load.assert_called_once_with()
        self.entrypoint2.load.assert_called_once_with()

        self.plugin1.assert_called_once_with()
        self.plugin2.assert_called_once_with()

        self.assertIn(self.plugin1.name, em.exts)
        self.assertIn(self.plugin2.name, em.exts)
        self.assertEqual(self.plugin1.return_value, em.exts['plugin1'])
        self.assertEqual(self.plugin2.return_value, em.exts['plugin2'])

    def test_configure_extensions_calls_configure_on_all_extensions(self):
        em = ExtensionManager(group='myplugins')
        em.load_extensions()
        em.configure_extensions('foo', bar='baz')

        self.plugin1.return_value.configure.assert_called_once_with(
            'foo', bar='baz'
        )

        self.plugin2.return_value.configure.assert_called_once_with(
            'foo', bar='baz'
        )

    def test_close_extensions_calls_close_on_all_extensions(self):
        em = ExtensionManager(group='myplugins')
        em.load_extensions()
        em.close_extensions()

        self.plugin1.return_value.close.assert_called_once_with()
        self.plugin2.return_value.close.assert_called_once_with()

    def test_setup_argument_parsers_calls_setup_argument_parser_on_all_extensions(self):
        em = ExtensionManager(group='myplugins')
        em.load_extensions()
        em.setup_argument_parsers(sentinel.parser, sentinel.subparsers)

        self.plugin1.return_value.setup_argument_parser.assert_called_once_with(
            sentinel.parser, sentinel.subparsers
        )
        self.plugin2.return_value.setup_argument_parser.assert_called_once_with(
            sentinel.parser, sentinel.subparsers
        )

    def test_loading_two_extensions_with_the_same_name_raises_exception(self):
        self.plugin1.name = 'fooplugin'
        self.plugin2.name = 'fooplugin'

        em = ExtensionManager(group='myplugins')
        with self.assertRaises(GWFError):
            em.load_extensions()
