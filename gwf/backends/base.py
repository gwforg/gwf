import abc
import json
import os.path
import logging
from collections import UserDict

from ..exceptions import NoLogFoundError

logger = logging.getLogger(__name__)


class Backend:

    """Abstract base class for backends."""

    def __init__(self, working_dir):
        self.working_dir = working_dir

    @property
    def supported_options(self):
        """Return the options supported on targets."""
        return set()

    @property
    def option_defaults(self):  # pragma: no cover
        """Return defaults for required target options."""
        return {}

    @abc.abstractmethod
    def submitted(self, target):
        """Return whether the target has been submitted."""

    @abc.abstractmethod
    def running(self, target):
        """Return whether the target is running."""

    @abc.abstractmethod
    def completed(self, target):
        """Return whether the target has completed."""

    @abc.abstractmethod
    def submit(self, target, dependencies):
        """Submit a target."""

    @abc.abstractmethod
    def cancel(self, target):
        """Cancel a target."""

    def logs(self, target, stderr=False):
        """Return log files for a target.

        If the backend cannot return logs a
        :class:`~gwf.exceptions.NoLogFoundError` is raised.

        By default standard output (stdout) is returned. If `stderr=True`
        standard error will be returned instead.

        :param gwf.Target target:
            Target to return logs for.
        :param bool stderr:
            default: False. If true, return standard error.
        :return:
            A file-like object. The user is responsible for closing the
            returned file(s) after use.
        :raises gwf.exceptions.NoLogFoundError:
            if the backend could not find a log for the given target.
        """
        try:
            if stderr:
                return open(self._log_path(target, 'stderr'), 'r')
            return open(self._log_path(target, 'stdout'), 'r')
        except OSError as e:
            raise NoLogFoundError()

    def _log_dir(self, target):
        """Path to directory containing logs for `target`."""
        return os.path.join(self.working_dir, '.gwf', 'logs')

    def _log_path(self, target, extension):
        return os.path.join(self._log_dir(target), '{}.{}'.format(target.name, extension))

    def close(self):
        """Close the backend."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class PersistableDict(UserDict):
    """A dictionary which can persist itself to JSON."""

    def __init__(self, path):
        super().__init__()

        self.path = path
        try:
            with open(self.path) as fileobj:
                self.data.update(json.load(fileobj))
        except (OSError, ValueError):
            # Catch ValueError for compatibility with Python 3.4.2. I haven't been
            # able to figure out what is different between 3.4.2 and 3.5 that
            # causes this. Essentially, 3.4.2 raises a ValueError saying that it
            # cannot parse the empty string instead of raising an OSError
            # (FileNotFoundError does not exist in 3.4.2) saying that the file does
            # not exist.
            pass

    def persist(self):
        with open(self.path + '.new', 'w') as fileobj:
            json.dump(self.data, fileobj)
            fileobj.flush()
            os.fsync(fileobj.fileno())
            fileobj.close()
        os.rename(self.path + '.new', self.path)
