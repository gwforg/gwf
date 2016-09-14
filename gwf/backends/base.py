import abc
import logging

from ..core import PreparedWorkflow
from ..exceptions import WorkflowNotPreparedError
from ..ext import Extension
from ..utils import dfs

logger = logging.getLogger(__name__)


class Backend(Extension):

    """Base class for backends.

    A backend is initialized with an instance of
    :class:`gwf.core.PreparedWorkflow`.

    .. warning::
      You should never override `__init__` in subclasses of :class:`Backend`
      unless you really know what you're doing.

    :cvar supported_options list: A list of the names of supported options.
    :cvar defaults dict: A dictionary with option defaults.
    """

    def __init__(self, workflow):
        if not isinstance(workflow, PreparedWorkflow):
            raise WorkflowNotPreparedError()
        self.workflow = workflow

        all_options = {option_name
                       for target in workflow.targets.values()
                       for option_name in target.options}

        for option_name in all_options:
            if option_name not in self.supported_options:
                logger.warn(
                    'Backend "%s" does not support option "%s".',
                    self.name,
                    option_name
                )

    @property
    @abc.abstractmethod
    def name(self):
        """User-friendly, single word name of the backend."""

    @property
    def supported_options(self):
        """Return the options supported on targets."""

    @property
    def option_defaults(self):
        """Return defaults for required target options."""
        return {}

    @abc.abstractmethod
    def configure(self, **options):
        """Configure the backend.

        This method *must* be called before any other method on the backend
        is used. Unless the backend is initialized directly, *gwf* is
        responsible for calling :func:`configure` to configure the backend.
        """

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

    def close(self):
        """Close the backend."""

    def schedule(self, target):
        """Schedule and submit a :class:`Target`s and its dependencies."""
        logger.debug('Scheduling target %s.', target.name)

        if self.submitted(target):
            logger.debug('Target %s has already been submitted.', target.name)
            return []

        scheduled = []
        for dependency in dfs(target, self.workflow.dependencies):
            logger.info(
                'Scheduling dependency %s of %s',
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

            logger.info('Submitting dependency %s', dependency.name)

            self.submit(dependency)
            scheduled.append(dependency)

        return scheduled

    def schedule_many(self, targets):
        """Schedule a list of :class:`Target`s and their dependencies.

        Will schedule the targets in `targets` with :func:`schedule`
        and return a list of schedules.

        :param list targets: A list of targets to be scheduled.
        :return: A list of schedules, one for each target in `targets`.
        """
        schedules = []
        for target in targets:
            schedules.append(self.schedule(target))
        return schedules
