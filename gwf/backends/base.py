from ..core import PreparedWorkflow
from ..exceptions import GWFError, WorkflowNotPreparedError
from ..utils import cache

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

    A workflow is initialized with a :class:`gwf.core.Workflow`
    """

    def __init__(self, workflow):
        if not isinstance(workflow, PreparedWorkflow):
            raise WorkflowNotPreparedError()
        self.workflow = workflow

    def submitted(self, targets):
        raise NotImplementedError()

    def running(self, targets):
        raise NotImplementedError()

    def submit(self, targets):
        raise NotImplementedError()

    def cancel(self, targets):
        raise NotImplementedError()

    def schedule(self, target):
        """Linearize the targets to be run.

        Returns a list of `Target`s to be run (in the order they should
        be submitted to the backend to make sure dependencies are handled
        correctly) and a set of the names of tasks that will be scheduled
        (to make sure dependency flags are set in the submission command).
        """

        root = target

        # If the target is already in the queue we just dismiss the scheduling
        # right away... this because we need to handle dependent nodes in the
        # queue differently, since for those we need to wait for completion.
        if self.submitted(root):
            return [], set()

        processed = set()
        scheduled = set()
        job_schedule = []

        def dfs(target):
            if target in processed:
                # we have already processed the node, and
                # if we should run the target name is scheduled
                # otherwise it isn't.
                return target in scheduled

            # schedule all dependencies before we schedule this task
            dependencies = self.workflow.dependencies[target]
            for dep in dependencies:
                dfs(dep)

            # If this task needs to run, then schedule it
            if self.submitted(target) or self.should_run(target):
                job_schedule.append(target)
                scheduled.add(target)

            processed.add(target)

        dfs(root)

        return job_schedule, scheduled

    def should_run(self, target):
        if any(self.should_run(dep) for dep in self.workflow.dependencies[target]):
            return True

        if target.is_sink:
            return True

        if any(self.workflow.file_cache[path] is None for path in target.outputs):
            return True

        if target.is_source:
            return False

        youngest_in_ts, youngest_in_path = max(self.workflow.file_cache[path]
                                               for path in target.inputs)
        oldest_out_ts, oldest_out_path = min(self.workflow.file_cache[path]
                                             for path in target.outputs)

        return youngest_in_ts > oldest_out_ts
