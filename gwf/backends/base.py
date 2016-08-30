import logging

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

        if not hasattr(cls, 'name'):
            raise GWFError(
                'Backend {} does not declare name class variable.'.format(name)
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

    def schedule(self, target):
        """Schedule and submit a :class:`Target`s and its dependencies."""
        if self.submitted(target):
            return []

        scheduled = []
        for scheduled_target in dfs(target, self.workflow.dependencies):
            logger.info('Scheduling target %s', scheduled_target.name)
            if not self.submitted(target) and self.workflow.should_run(target):
                logger.info('Submitting target %s', scheduled_target.name)
                self.submit(scheduled_target)
                scheduled.append(scheduled_target)
        return scheduled
