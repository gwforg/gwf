from gwf.exceptions import GWFError


class BackendError(GWFError):
    """Base class for backend errors."""


class UnsupportedOperationError(BackendError):
    """Operation not supported by this backend."""


class UnknownDependencyError(BackendError):
    pass


class UnknownTargetError(BackendError):
    pass


class LogNotFoundError(BackendError):
    """No log found for target."""

    def __init__(self):
        super().__init__('No log found.')