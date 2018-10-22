from gwf.exceptions import GWFError


class BackendError(GWFError):
    pass


class UnsupportedOperationError(BackendError):
    pass


class DependencyError(BackendError):
    pass


class TargetError(BackendError):
    pass


class LogError(BackendError):
    def __init__(self):
        super().__init__("Log not found.")
