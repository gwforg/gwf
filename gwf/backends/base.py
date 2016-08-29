from ..core import PreparedWorkflow
from ..exceptions import GWFError, WorkflowNotPreparedError
from ..utils import dfs

BACKENDS = {}


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

    def submitted(self, targets):
        """Return whether the target has been submitted."""
        raise NotImplementedError()

    def running(self, targets):
        """Return whether the target is running."""
        raise NotImplementedError()

    def submit(self, targets):
        """Submit a target."""
        raise NotImplementedError()

    def cancel(self, targets):
        """Cancel a target."""
        raise NotImplementedError()

    def _get_schedule_for_target(self, root):
        """Linearize the targets to be run.

        Returns a list of :class:`~gwf.Target`s to be run (in the order they
        should be submitted to the backend to make sure dependencies are handled
        correctly) and a set of the names of tasks that will be scheduled
        (to make sure dependency flags are set in the submission command).
        """

        # If the target is already in the queue we just dismiss the scheduling
        # right away... this because we need to handle dependent nodes in the
        # queue differently, since for those we need to wait for completion.
        if self.submitted(root):
            return [], set()

        job_schedule = []
        for target in dfs(root, self.workflow.dependencies):
            if self.submitted(target) or self.workflow.should_run(target):
                job_schedule.append(target)
        return job_schedule, set(job_schedule)
