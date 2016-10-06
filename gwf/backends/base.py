import abc
import logging

from ..events import post_schedule, pre_schedule
from ..ext import Extension
from ..utils import dfs

logger = logging.getLogger(__name__)


class Backend(Extension):

    """Abstract base class for backends."""

    @property
    def supported_options(self):
        """Return the options supported on targets."""

    @property
    def option_defaults(self):  # pragma: no cover
        """Return defaults for required target options."""
        return {}

    def configure(self, workflow, config):
        """Configure the backend.

        This method *must* be called before any other method on the backend
        is used. Unless the backend is initialized directly, *gwf* is
        responsible for calling :func:`configure` to configure the backend.
        """
        self.workflow = workflow
        self.config = config

    @abc.abstractmethod
    def submitted(self, target):
        """Return whether the target has been submitted."""

    @abc.abstractmethod
    def running(self, target):
        """Return whether the target is running."""

    @abc.abstractmethod
    def submit(self, target):
        """Submit a target."""

    @abc.abstractmethod
    def cancel(self, target):
        """Cancel a target."""

    @abc.abstractmethod
    def logs(self, target, stderr=False):
        """Return log files for a target.

        If the backend cannot return logs a
        :class:`~gwf.exceptions.NoLogFoundError` is raised.

        By default only standard output (stdout) is returned. If `stderr=True`
        the function will return a tuple (stdout, stderr).

        :param gwf.Target target:
            Target to return logs for.
        :param bool stderr:
            default: False. If true, the method will return a tuple consisting
            of both the standard and error output.
        :return:
            A file-like object or a tuple (stdout, stderr) of file-like objects.
            The user is responsible for closing the returned file(s) after use.
        :raises gwf.exceptions.NoLogFoundError:
            if the backend could not find a log for the given target.
        """

    def close(self):
        """Close the backend."""

    def schedule(self, target):
        """Schedule and submit a :class:`gwf.Target` and its dependencies.

        This method is provided by :class:`Backend` and should not be overriden.
        """
        logger.info('Scheduling target %s.', target.name)

        if self.submitted(target):
            logger.debug('Target %s has already been submitted.', target.name)
            return []

        scheduled = []
        for dependency in dfs(target, self.workflow.dependencies):
            if dependency.name != target.name:
                logger.info(
                    'Scheduling dependency %s of %s.',
                    dependency.name,
                    target.name
                )

            if self.submitted(dependency):
                logger.debug(
                    'Target %s has already been submitted.',
                    dependency.name
                )
                continue

            if not self.workflow.should_run(dependency):
                logger.debug(
                    'Target %s should not run.',
                    dependency.name
                )
                continue

            logger.info('Submitting target %s.', dependency.name)

            self.submit(dependency)
            scheduled.append(dependency)

        return scheduled

    def schedule_many(self, targets):
        """Schedule a list of :class:`gwf.Target` and their dependencies.

        Will schedule the targets in `targets` with :func:`schedule`
        and return a list of schedules.

        This method is provided by :class:`Backend` and should not be overriden.

        :param list targets: A list of targets to be scheduled.
        :return: A list of schedules, one for each target in `targets`.
        """
        pre_schedule.trigger(targets=targets)

        schedules = []
        for target in targets:
            schedules.append(self.schedule(target))

        post_schedule.trigger(targets=targets, schedules=schedules)
        return schedules
