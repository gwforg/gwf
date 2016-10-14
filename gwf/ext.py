import abc
import logging

from pkg_resources import DistributionNotFound, iter_entry_points

from .exceptions import GWFError

logger = logging.getLogger(__name__)


class Extension(abc.ABC):

    def configure(self, *args, **kwargs):
        """Called to configure the extension."""

    def close(self):
        """Called when the extension is closed."""

    def setup_argument_parser(self, parser, subparsers):
        """Modify the main argument parser and subparsers.

        This method is called with an instance of
        :class:`argparse.ArgumentParser` and the extension may add any
        subcommands and arguments to the parser as long as they don't conflict
        with other subcommands/arguments.
        """

    def setup_subparser(self, subparsers, name, description, handler):
        """Helper method for setting up subparsers.

        :param subparsers:
            A subparsers object as returned by
            :func:`ArgumentParser.add_subparsers`.
        :param name:
            Name of the subcommand.
        :param description:
            A short description of the subcommand.
        :param handler:
            The callable that will be called when the subcommand is used.

        This method is most often used in :func:`setup_argument_parser` to set
        up subcommands. For example, to register a subcommand called ``foo``
        that calls the method ``on_foo``::

            class MyPlugin(Plugin):
                name = 'myplugin'

                def setup_argument_parser(self, parser, subparsers):
                    self.setup_subparser(
                        subparsers,
                        'foo',
                        'a command that prints a greeting.',
                        self.on_foo,
                    )

                def on_foo(self):
                    print('Hello!')

        When the plugin is installed and registered properly, the user will now
        be able to run ``gwf foo`` which will print ``Hello!`` and exit.
        """
        subparser = subparsers.add_parser(name, description=description)
        subparser.set_defaults(func=handler)
        return subparser


class ExtensionManager:

    def __init__(self, group):
        self.group = group
        self.exts = {}

    def load_extensions(self):
        for entry_point in iter_entry_points(group=self.group, name=None):
            ext_cls = entry_point.load()
            if entry_point.name in self.exts:
                raise GWFError(
                    'Extension with name "{}" already loaded.'.format(
                        entry_point.name
                    )
                )
            self.exts[entry_point.name] = ext_cls()

    def configure_extensions(self, *args, **kwargs):
        for ext in self.exts.values():
            ext.configure(*args, **kwargs)

    def close_extensions(self):
        for ext in self.exts.values():
            ext.close()

    def setup_argument_parsers(self, parser, subparsers):
        for ext in self.exts.values():
            ext.setup_argument_parser(parser, subparsers)
