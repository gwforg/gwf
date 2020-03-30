import logging
from enum import Enum
from pkg_resources import iter_entry_points

from ..utils import PersistableDict, retry
from .exceptions import BackendError, DependencyError, TargetError
from .logmanager import FileLogManager

logger = logging.getLogger(__name__)


__all__ = ("Backend", "Status")


def _load_backends():
    return {ep.name: ep.load() for ep in iter_entry_points("gwf.backends")}


class Status(Enum):
    """Status of a target.

    A target is unknown to the backend if it has not been submitted or the
    target has completed and thus isn't being tracked anymore by the backend.

    A target is submitted if it has been successfully submitted to the backend
    and is pending execution.

    A target is running if it is currently being executed by the backend.
    """

    UNKNOWN = 0  #: The backend is not aware of the status of this target (it may be completed or failed).
    SUBMITTED = 1  #: The target has been submitted, but is not currently running.
    RUNNING = 2  #: The target is currently running.


class Backend:
    """Base class for backends."""

    option_defaults = {}
    log_manager = FileLogManager()

    @classmethod
    def list(cls):
        """Return the names of all registered backends."""
        return set(_load_backends().keys())

    @classmethod
    def from_name(cls, name):
        """Return backend class for the backend given by `name`.

        Returns the backend class registered with `name`. Note that the *class*
        is returned, not the instance, since not all uses requires
        initialization of the backend (e.g. accessing the backends' log
        manager), and initialization of the backend may be expensive.

        :arg str name: Path to a workflow file, optionally specifying a
            workflow object in that file.
        """
        return _load_backends()[name]

    @classmethod
    def from_config(cls, config):
        """Return backend class for the backend specified by `config`.

        See :func:`Backend.from_name` for further information."""
        return cls.from_name(config["backend"])

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def status(self, target):
        """Return the status of `target`.

        :param gwf.Target target: The target to return the status of.
        :return gwf.backends.Status: Status of `target`.
        """

    def submit(self, target, dependencies):
        """Submit `target` with `dependencies`.

        This method must submit the `target` and return immediately. That is,
        the method must not block while waiting for the target to complete.

        :param gwf.Target target:
            The target to submit.
        :param dependencies:
            An iterable of :class:`gwf.Target` objects that `target` depends on
            and that have already been submitted to the backend.
        """

    def cancel(self, target):
        """Cancel `target`.

        :param gwf.Target target:
            The target to cancel.
        :raises gwf.exception.TargetError:
            If the target does not exist in the workflow.
        """

    @classmethod
    def logs(cls, target, stderr=False):
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
        if stderr:
            return cls.log_manager.open_stderr(target)
        return cls.log_manager.open_stdout(target)

    def close(self):
        """Close the backend.

        Called when the backend is no longer needed and should close all
        resources (open files, connections) used by the backend.
        """


class PbsLikeBackendBase(Backend):
    """PBS-like backend base class."""

    option_flags = {}
    option_defaults = {}
    log_manager = FileLogManager()

    def __init__(self):
        try:
            self._status = self.parse_queue_output(self.call_queue_command())
        except retry.RetryError as exc:
            raise BackendError("Could not get queue state") from exc

        class_name = self.__class__.__name__
        backend_name = class_name.strip("Backend").lower()

        path = ".gwf/{name}-backend-tracked.json".format(name=backend_name)
        self._tracked = PersistableDict(path=path)

    def parse_queue_output(self):
        raise NotImplementedError("parse_queue_output")

    def call_queue_command(self):
        raise NotImplementedError("call_queue_command")

    def call_submit_command(self):
        raise NotImplementedError("call_submit_command")

    def call_cancel_command(self):
        raise NotImplementedError("call_cancel_command")

    def compile_script(self, target):
        raise NotImplementedError("compile_script")

    def status(self, target):
        try:
            return self._get_status(target)
        except KeyError:
            return Status.UNKNOWN

    def submit(self, target, dependencies):
        script = self.compile_script(target)
        dependency_ids = self._collect_dependency_ids(dependencies)
        try:
            stdout = self.call_submit_command(script, dependency_ids)
        except retry.RetryError as exc:
            raise BackendError("Could not submit target") from exc
        else:
            job_id = stdout.strip()
            self._add_job(target, job_id)

    def cancel(self, target):
        try:
            job_id = self.get_job_id(target)
            self.call_cancel_command(job_id)
        except KeyError as exc:
            raise TargetError(target.name) from exc
        except retry.RetryError as exc:
            raise BackendError("Could not cancel target") from exc
        else:
            self.forget_job(target)

    def close(self):
        self._tracked.persist()

    def forget_job(self, target):
        """Force the backend to forget the job associated with `target`."""
        job_id = self.get_job_id(target)
        del self._status[job_id]
        del self._tracked[target.name]

    def get_job_id(self, target):
        """Get the Slurm job id for a target.

        :raises KeyError: if the target is not tracked by the backend.
        """
        return self._tracked[target.name]

    def _add_job(self, target, job_id, initial_status=Status.SUBMITTED):
        self._set_job_id(target, job_id)
        self._set_status(target, initial_status)

    def _set_job_id(self, target, job_id):
        self._tracked[target.name] = job_id

    def _get_status(self, target):
        job_id = self.get_job_id(target)
        return self._status[job_id]

    def _set_status(self, target, status):
        job_id = self.get_job_id(target)
        self._status[job_id] = status

    def _collect_dependency_ids(self, dependencies):
        try:
            return [self._tracked[dep.name] for dep in dependencies]
        except KeyError as exc:
            raise DependencyError(exc.args[0])
