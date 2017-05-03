import json
import os.path
import logging
import warnings
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

    UNKNOWN = 0
    SUBMITTED = 1
    RUNNING = 2


class UnknownDependencyError(BackendError):
    pass


class UnknownTargetError(BackendError):
    pass


def inherit_options(func, super_options):
    @wraps(func)
    def inner_wrapper(self, target, *args, **kwargs):
        target.inherit_options(super_options)
        return func(self, target, *args, **kwargs)
    return inner_wrapper


def check_options(func, supported_options):
    @wraps(func)
    def inner_wrapper(self, target, *args, **kwargs):
        for option_name in list(target.options.keys()):
            if option_name not in supported_options:
                warnings.warn(
                    'Backend does not support option {} used in {}. Option will be ignored.'.format(
                        option_name, target.name
                    )
                )
                del target.options[option_name]
        return func(self, target, *args, **kwargs)
    return inner_wrapper


class BackendType(type):
    """A metaclass for backends.

    All backends are initialized via this metaclass.
    """
    def __new__(metacls, name, bases, namespace, **kwargs):
        # Check that all required methods exist. The logs() method isn't required,
        # since a default implementation is provided by Backend.
        for method_name in ('submit', 'cancel', 'status', 'close'):
            if method_name not in namespace:
                raise BackendError('Invalid backend implementation. Backend does not implement {}.'.format(method_name))

        if not bases:
            return type.__new__(metacls, name, bases, namespace)

        # Decorate the submit method with a decorator that injects backend
        # target defaults into the target options.
        option_defaults = namespace.get('option_defaults', {})
        namespace['submit'] = inherit_options(namespace['submit'], option_defaults)

        # Decorate the submit method with a decorator that checks whether the
        # option in the given target are supported by the backend. Warns the
        # user and removes the option if this is not the case.
        supported_options = namespace.get('supported_options', {})
        namespace['submit'] = check_options(namespace['submit'], supported_options)

        return type.__new__(metacls, name, bases, namespace)


class Backend(metaclass=BackendType):
    """Abstract base class for backends."""

    def __init__(self, working_dir):
        self.working_dir = working_dir

    def status(self, target):
        """Return the status of a target."""

    def submit(self, target, dependencies):
        """Submit a target."""

    def cancel(self, target):
        """Cancel a target."""

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

    def _log_dir(self, target):
        """Path to directory containing logs for `target`."""
        return os.path.join(self.working_dir, '.gwf', 'logs')

    def _log_path(self, target, extension):
        return os.path.join(self._log_dir(target), '{}.{}'.format(target.name, extension))

    def close(self):
        """Close the backend."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


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
