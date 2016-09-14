import abc
import logging

from pkg_resources import iter_entry_points

logger = logging.getLogger(__name__)


class Extension(abc.ABC):

    @property
    @abc.abstractmethod
    def name(self):
        pass

    @abc.abstractmethod
    def configure(self, *args, **kwargs):
        pass

    @abc.abstractmethod
    def close(self, *args, **kwargs):
        pass

    def setup_argument_parser(self, parser):
        """Modify the main argument parser.

        This static method is called with an instance of
        :class:`argparse.ArgumentParser` and the backend may add any
        subcommands and arguments to the parser as long as they don't conflict
        with other subcommands/arguments.
        """
        pass


class Plugin(Extension):
    pass


class ExtensionManager:

    def __init__(self):
        self.extensions = self.reload()

    def reload(self):
        return [
            entry_point.load()
            for entry_point in iter_entry_points(group='gwf.ext', name=None)
        ]

    def list(self, parent=None):
        if parent is None:
            return list(self.extensions)
        return [ext for ext in self.extensions if issubclass(ext, parent)]
