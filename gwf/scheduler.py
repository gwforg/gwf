import os.path
from collections import defaultdict

from .exceptions import GWFError
from .utils import cache, iter_inputs, iter_outputs


_ex_msg_file_provided_by_multiple_targets = (
    'File "{}" provided by multiple targets "{}" and "{}".'
)

_ex_msg_file_required_but_not_provided = (
    'File "{}" is required by "{}", but does not exist and is not provided by '
    'a target.'
)

_ex_msg_target_depends_on_itself = (
    'Target {} depends on itself.'
)


@cache
def _get_file_timestamp(filename):
    return os.path.getmtime(filename)


def _compute_dependency_graph(targets):
    provides = {}
    dependencies = defaultdict(list)
    dependents = defaultdict(list)

    for target, path in iter_outputs(targets):
        if path in provides:
            raise GWFError(
                _ex_msg_file_provided_by_multiple_targets.format(
                    path, provides[path].name, target.name
                )
            )

        provides[path] = target

    for target, path in iter_inputs(targets):
        if os.path.exists(path):
            continue

        if path not in provides:
            raise GWFError(
                _ex_msg_file_required_but_not_provided.format(
                    path, target.name
                )
            )

        dependencies[target].append(provides[path])

    for target, deps in dependencies.items():
        for dep in deps:
            dependents[dep].append(target)

    return provides, dependencies, dependents


class Scheduler:
    """Run a workflow on a specific backend.

    Three attributes are available on the object:

    `provides`:
        a dictionary of file paths and the corresponding targets that
        produce the files.

    `dependencies`:
        a dictionary where each key is a `Target` and the value is a list
        of `Target`s that the key `Target` depends on.

    `dependents`:
        a dictionary where each key is a `Target` and the value is a list
        of `Targets`s that depend on the key `Target`.
    """

    def __init__(self, workflow, backend):
        self.workflow = workflow
        self.backend = backend

        self.provides, self.dependencies, self.dependents = \
            _compute_dependency_graph(self.workflow.targets.values())

        self._check_circular_dependencies()

    def _check_circular_dependencies(self):
        for target in self.workflow.targets.values():
            for dep in self.dependencies[target]:
                if target in self._get_dependencies(dep):
                    raise GWFError(
                        _ex_msg_target_depends_on_itself.format(target.name)
                    )

    @cache
    def should_run(self, target):
        if any(self.should_run(dep) for dep in self.dependencies[target]):
            return True

        if target.is_sink:
            return True

        if any(not os.path.exists(path) for path in target.outputs):
            return True

        if target.is_source:
            return False

        youngest_in_ts, youngest_in_path = \
            max((_get_file_timestamp(path), path) for path in target.inputs)

        oldest_out_ts, oldest_out_path = \
            min((_get_file_timestamp(path), path) for path in target.outputs)

        return youngest_in_ts > oldest_out_ts

    def _get_dependencies(self, target):
        """Return all tasks necessary for building the target.

        The set of tasks is just returned as set.
        """
        # Working with a list to preserve the order. It makes lookups slower
        # but hopefully these sets won't be terribly long ... if it becomes a
        # problem it is easy enough to fix it.
        processed = []

        def dfs(other_target):
            if other_target in processed:
                return
            else:
                processed.append(other_target)
                for dep in self.dependencies[other_target]:
                    dfs(dep)

        dfs(target)
        return processed

    def _get_execution_schedule_for_target(self, target):
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
        if self.backend.submitted(root):
            return [], set()

        processed = set()
        scheduled = set()
        job_schedule = []

        def dfs(target):
            if target in processed:
                # we have already processed the node, and
                # if we should run the target name is scheduled
                # otherwise it isn't.
                return target.name in scheduled

            # schedule all dependencies before we schedule this task

            # TODO: check whether it is possible to just use
            #       dependencies[target] here.
            for dep in self.dependencies[target]:
                dfs(dep)

            # If this task needs to run, then schedule it
            if self.backend.submitted(target) or self.should_run(target):
                job_schedule.append(target)
                scheduled.add(target)

            processed.add(target)

        dfs(root)

        return job_schedule, scheduled

    def run(self):
        for target in self.workflow.targets.values():
            schedule, scheduled = self._get_execution_schedule_for_target
