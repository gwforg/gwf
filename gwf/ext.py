import abc
import logging

from pkg_resources import iter_entry_points

logger = logging.getLogger(__name__)


class Extension(abc.ABC):

    @property
    @abc.abstractmethod
    def name(self):
        """Single-word, user-friendly name for the extension."""

    @abc.abstractmethod
    def configure(self, *args, **kwargs):
        """Called to configure the extension."""

    @abc.abstractmethod
    def close(self, *args, **kwargs):
        """Called when the extension is closed."""

    @classmethod
    def setup_argument_parser(cls, parser):
        """Modify the main argument parser.

        This static method is called with an instance of
        :class:`argparse.ArgumentParser` and the extension may add any
        subcommands and arguments to the parser as long as they don't conflict
        with other subcommands/arguments.
        """


class Plugin(Extension):
    pass


class ExtensionManager:

    group = 'gwf.ext'

    def __init__(self):
        self.extensions = self.reload()

    def reload(self):
        return [
            entry_point.load()
            for entry_point in iter_entry_points(group=self.group, name=None)
        ]

    def list(self, parent=None):
        if parent is None:
            return list(self.extensions)
        return [ext for ext in self.extensions if issubclass(ext, parent)]
