import json
import logging
import os.path
from enum import Enum

import attrs

from ..utils import entry_points
from .exceptions import TargetError

logger = logging.getLogger(__name__)


class BackendStatus(Enum):
    """BackendStatus of a target.

    A target is unknown to the backend if it has not been submitted or the
    target has completed and thus isn't being tracked anymore by the backend.

    A target is submitted if it has been successfully submitted to the backend
    and is pending execution.

    A target is running if it is currently being executed by the backend.
    """

    UNKNOWN = 0
    SUBMITTED = 1
    RUNNING = 2
    COMPLETED = 3
    FAILED = 4


def guess_backend():
    max_score = -1000
    chosen_backend = None
    for backend_name, (_, score) in discover_backends().items():
        if score > max_score:
            max_score = score
            chosen_backend = backend_name
    return max_score, chosen_backend


def discover_backends():
    return {ep.name: ep.load() for ep in entry_points(group="gwf.backends")}


def list_backends():
    """Return the names of all registered backends."""
    return set(discover_backends().keys())


def create_backend(name, working_dir, config):
    """Return backend class for the backend given by `name`.

    Returns the backend class registered with `name`. Note that the *class*
    is returned, not the instance, since not all uses requires
    initialization of the backend (e.g. accessing the backends' log
    manager), and initialization of the backend may be expensive.

    :arg str name: Path to a workflow file, optionally specifying a
        workflow object in that file.
    """
    backend_args = config.get_namespace(f"backend.{name}")
    backend_cls, _ = discover_backends()[name]
    return backend_cls(working_dir=working_dir, **backend_args)


@attrs.define()
class TrackingBackend:
    working_dir: str = attrs.field()
    name: str = attrs.field()
    ops: object = attrs.field()

    _tracked_jobs: dict = attrs.field(init=False, repr=False)
    _job_states: dict = attrs.field(init=False, repr=False)

    @_tracked_jobs.default
    def _init_tracked(self):
        try:
            with open(self._get_state_path()) as state_file:
                return json.load(state_file)
        except FileNotFoundError:
            return {}

    @_job_states.default
    def _init_status(self):
        return self.ops.get_job_states(self._tracked_jobs.values())

    def _get_state_path(self):
        return os.path.join(
            self.working_dir, ".gwf", f"{self.name}-backend-tracked.json"
        )

    def status(self, target):
        job_id = self._tracked_jobs.get(target.name)
        return self._job_states.get(job_id, BackendStatus.UNKNOWN)

    def submit(self, target, dependencies):
        dependency_ids = [self._tracked_jobs[dep.name] for dep in dependencies]
        job_id = self.ops.submit_target(target, dependency_ids)
        self._tracked_jobs[target.name] = job_id
        self._job_states[job_id] = BackendStatus.SUBMITTED

    def cancel(self, target):
        try:
            job_id = self._tracked_jobs[target.name]
            self.ops.cancel_job(job_id)
            del self._job_states[job_id]
            del self._tracked_jobs[target.name]
        except KeyError as exc:
            raise TargetError(target.name) from exc

    def close(self):
        self.ops.close()
        with open(self._get_state_path(), "w") as state_file:
            json.dump(self._tracked_jobs, state_file)

    @property
    def target_defaults(self):
        return self.ops.target_defaults

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
