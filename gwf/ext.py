import abc
import logging

from pkg_resources import iter_entry_points

logger = logging.getLogger(__name__)


class Extension(abc.ABC):

    @property
    @abc.abstractmethod
    def name(self):
        """Single-word, user-friendly name for the extension."""

    def configure(self, *args, **kwargs):
        """Called to configure the extension."""

    def close(self):
        """Called when the extension is closed."""

    def setup_argument_parser(self, parser, subparsers):
        """Modify the main argument parser and subparsers.

        This static method is called with an instance of
        :class:`argparse.ArgumentParser` and the extension may add any
        subcommands and arguments to the parser as long as they don't conflict
        with other subcommands/arguments.
        """

    def setup_subparser(self, subparsers, name, description, handler):
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

            logger.debug(
                'Loading extension: %s (%s).',
                ext_cls.name,
                self.group
            )

            logger.debug(
                'Initializing extension: %s (%s).',
                ext_cls.name,
                self.group
            )
            self.exts[ext_cls.name] = ext_cls()

    def configure_extensions(self, *args, **kwargs):
        for ext in self.exts.values():
            ext.configure(*args, **kwargs)

    def close_extensions(self):
        for ext in self.exts.values():
            ext.close()

    def setup_argument_parsers(self, parser, subparsers):
        for ext in self.exts.values():
            ext.setup_argument_parser(parser, subparsers)
