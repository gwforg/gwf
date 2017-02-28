import abc
import logging

from ..ext import Extension

logger = logging.getLogger(__name__)


class Backend(Extension):

    """Abstract base class for backends."""

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
    def failed(self, target):
        """Return whether the target has failed."""

    @abc.abstractmethod
    def completed(self, target):
        """Return whether the target has completed."""

    @abc.abstractmethod
    def submit(self, target, dependencies):
        """Submit a target."""

    @abc.abstractmethod
    def cancel(self, target):
        """Cancel a target."""

    @abc.abstractmethod
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

    def close(self):
        """Close the backend."""
