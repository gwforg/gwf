import logging
from enum import Enum
from functools import wraps

from pkg_resources import iter_entry_points

from .exceptions import BackendError
from .logmanager import FileLogManager

logger = logging.getLogger(__name__)


__all__ = (
    'list_backends',
    'backend_from_name',
    'backend_from_config',
    'Backend',
    'Status',
)


def _load_backends():
    return {ep.name: ep.load() for ep in iter_entry_points('gwf.backends')}


def list_backends():
    """Return the names of all registered backends."""
    return set(_load_backends().keys())


def backend_from_name(name):
    """Return backend class for the backend given by `name`.

    Returns the backend class registered with `name`. Note that the *class* is returned, not the instance, since not all
    uses requires initialization of the backend (e.g. accessing the backends' log manager), and initialization of the
    backend may be expensive.

    :arg str name:
        Path to a workflow file, optionally specifying a workflow object in that file.
    """
    return _load_backends()[name]


def backend_from_config(config):
    """Return backend class for the backend specified by `config`.

    See :func:`backend_from_name` for further information."""
    return backend_from_name(config['backend'])


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

    log_manager = FileLogManager()

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

        Called when the backend is no longer needed and should close
        all resources (open files, connections) used by the backend.
        """
