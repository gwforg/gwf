import click


class GWFError(click.ClickException):
    """Base exception for all gwf exceptions."""


class InvalidNameError(GWFError):
    """Invalid name."""


class InvalidTypeError(GWFError):
    """Raised when an argument to a target has the wrong type."""


class TargetExistsError(GWFError):
    """Target already exists in the workflow."""

    def __init__(self, target):
        self.target = target

        message = (
            'Target "{}" already exists in workflow.'
        ).format(target.name)

        super(TargetExistsError, self).__init__(message)


class TargetNotFoundError(GWFError):
    """Target does not exist in the workflow."""

    def __init__(self, target_name):
        self.target_name = target_name

        message = (
            'Target "{}" is not found in the workflow.'
        ).format(target_name)

        super(TargetNotFoundError, self).__init__(message)


class MultipleProvidersError(GWFError):
    """File is provided by multiple targets.

    :param str path: path of the file.
    :param str gwf.Target: one of the targets providing the file.
    :param str gwf.Target: the other target providing the file.
    """

    def __init__(self, path, target, other_target):
        self.path = path
        self.target = target
        self.other_target = other_target

        message = (
            'File "{}" provided by targets "{}" and "{}".'
        ).format(self.path, self.target, self.other_target)

        super(MultipleProvidersError, self).__init__(message)


class MissingProviderError(GWFError):
    """File required by a target does not exist and is not provided by a target.

    :param str path:
        path to the file.
    :param gwf.Target target:
        target that requires the file.
    """

    def __init__(self, path, target):
        self.path = path
        self.target = target

        message = (
            'File "{}" is required by "{}", but does not exist and is not '
            'provided by any target in the workflow.'
        ).format(self.path, self.target)

        super(MissingProviderError, self).__init__(message)


class CircularDependencyError(GWFError):
    """Circular dependency.

    :param gwf.Target target:
        target that depends on itself.
    """

    def __init__(self, target):
        self.target = target

        message = (
            'Target {} depends on itself.'
        ).format(self.target)

        super(CircularDependencyError, self).__init__(message)


class IncludeWorkflowError(GWFError):
    """Workflow could not be included."""


class InvalidPathError(GWFError):
    """Target declared a directory as as input or output."""
