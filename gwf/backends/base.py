import logging
import warnings

from ..core import PreparedWorkflow
from ..exceptions import GWFError, WorkflowNotPreparedError
from ..utils import dfs

BACKENDS = {}

logger = logging.getLogger(__name__)


def register_backend(name, backend_cls):
    BACKENDS[name] = backend_cls


def get_backends():
    return dict(BACKENDS)


class BackendType(type):

    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)

        # Do not register the base Backend class.
        if name == 'Backend':
            return cls

        if 'name' not in class_dict:
            raise GWFError(
                'Backend {} does not declare name class variable.'.format(name)
            )

        if 'schedule' in class_dict or 'schedule_all' in class_dict:
            warnings.warn(
                'Subclasses of Backend should not override schedule() or '
                'schedule_all().'
            )

        register_backend(getattr(cls, 'name'), cls)
        return cls


class Backend(metaclass=BackendType):

    """Representation of a backend.

    A backend is initialized with an instance of
    :class:`gwf.core.PreparedWorkflow`.
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


    def configure(self, **options):
        """Configure the backend.

        This method *must* be called before any other method on the backend
        is used. Unless the backend is initialized directly, *gwf* is
        responsible for calling :func:`configure` to configure the backend.
        """
        pass

    def submitted(self, target):
        """Return whether the target has been submitted."""
        pass

    def running(self, target):
        """Return whether the target is running."""
        pass

    def submit(self, target):
        """Submit a target."""
        pass

    def cancel(self, target):
        """Cancel a target."""
        pass

    def close(self):
        """Close the backend."""
        pass

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
