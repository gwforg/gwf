import abc
import logging

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


class Plugin(Extension):
    pass
