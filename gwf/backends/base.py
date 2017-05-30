import json
import os.path
import logging
from collections import UserDict
from enum import Enum
from functools import wraps

from ..exceptions import NoLogFoundError, BackendError

logger = logging.getLogger(__name__)


class Status(Enum):
    """Status of a target.

    A target is unknown to the backend if it has not been submitted
    or the target has completed and thus isn't being tracked anymore by
    the backend.

    A target is submitted if it has been successfully submitted to
    the backend and is pending execution.

    A target is running if it is currently being executed by the backend.
    """

    UNKNOWN = 0  #: The backend is not aware of the status of this target (it may be completed or failed).
    SUBMITTED = 1  #: The target has been submitted, but is not currently running.
    RUNNING = 2  #: The target is currently running.


class UnknownDependencyError(BackendError):
    pass


class UnknownTargetError(BackendError):
    pass


def check_options(func, supported_options, super_options):
    """Decorator for the submit method of backends.

    Injects backend target defaults into the target options and
    checks whether the option in the given target are supported by
    the backend. Warns the user and removes the option if this is
    not the case.
    """
    @wraps(func)
    def inner_wrapper(self, target, *args, **kwargs):
        target.inherit_options(super_options)
        for option_name, option_value in list(target.options.items()):
            if option_name not in supported_options:
                logger.warning(
                    'Option "{}" used in "{}" is not supported by backend. Ignored.'.format(
                        option_name, target.name
                    )
                )
                del target.options[option_name]
            elif option_value is None:
                del target.options[option_name]
        return func(self, target, *args, **kwargs)
    return inner_wrapper


class BackendType(type):
    """A metaclass for backends.

    All backends are initialized via this metaclass. It has three distinct
    responsibilities:

        1. Check whether implementations of the Backend base class provide all
           required methods: submit, cancel, status and close.

        2. Wrap the submit method in subclasses with a decorator which modifies
           the given target's options with the option defaults provided by the
           backend through the ``option_defaults`` attribute. If the attribute
           is not supplied, submit will not be decorated.

        3. Wrap the submit method in subclasses with a decorator which removes
           options used in the given target that are not supported by the
           backend.

    This is not necessary magic, however, it removes a lot of boilerplate code
    from backend implementations and ensures consistency in warnings.
    """
    def __new__(mcs, name, bases, namespace, **kwargs):
        if not bases:
            return type.__new__(mcs, name, bases, namespace)

        # Check that all required methods exist. The logs() method isn't required,
        # since a default implementation is provided by Backend.
        for method_name in ('submit', 'cancel', 'status', 'close'):
            if method_name not in namespace:
                raise BackendError('Invalid backend implementation. Backend does not implement {}.'.format(method_name))

        option_defaults = namespace.get('option_defaults', {})
        namespace['submit'] = check_options(namespace['submit'], option_defaults.keys(), option_defaults)
        return type.__new__(mcs, name, bases, namespace)


class Backend(metaclass=BackendType):
    """Base class for backends."""

    def status(self, target):
        """Return the status of `target`.

        :param gwf.Target target: The target to return the status of.
        :return gwf.backends.Status: Status of `target`.
        """

    def submit(self, target, dependencies):
        """Submit `target` with `dependencies`.

        This method must submit the `target` and return immediately. That is, the method
        must not block while waiting for the target to complete.

        :param gwf.Target target:
            The target to submit.
        :param dependencies:
            An iterable of :class:`gwf.Target` objects that `target` depends on and that have
            already been submitted to the backend.
        """

    def cancel(self, target):
        """Cancel `target`.

        :param gwf.Target target:
            The target to cancel.
        :raises gwf.exception.UnknownTargetError:
            If the target does not exist in the workflow.
        """

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
        try:
            if stderr:
                return open(self._log_path(target, 'stderr'), 'r')
            return open(self._log_path(target, 'stdout'), 'r')
        except OSError:
            raise NoLogFoundError()

    def _log_dir(self):
        """Path to directory containing logs for `target`."""
        return os.path.join('.gwf', 'logs')

    def _log_path(self, target, extension):
        return os.path.join(self._log_dir(), '{}.{}'.format(target.name, extension))

    def close(self):
        """Close the backend.

        Called when the backend is no longer needed and should close
        all resources (open files, connections) used by the backend.
        """


class PersistableDict(UserDict):
    """A dictionary which can persist itself to JSON."""

    def __init__(self, path):
        super().__init__()

        self.path = path
        try:
            with open(self.path) as fileobj:
                self.data.update(json.load(fileobj))
        except (OSError, ValueError):
            # Catch ValueError for compatibility with Python 3.4.2. I haven't been
            # able to figure out what is different between 3.4.2 and 3.5 that
            # causes this. Essentially, 3.4.2 raises a ValueError saying that it
            # cannot parse the empty string instead of raising an OSError
            # (FileNotFoundError does not exist in 3.4.2) saying that the file does
            # not exist.
            pass

    def persist(self):
        with open(self.path + '.new', 'w') as fileobj:
            json.dump(self.data, fileobj)
            fileobj.flush()
            os.fsync(fileobj.fileno())
            fileobj.close()
        os.rename(self.path + '.new', self.path)
