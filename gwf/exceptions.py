from __future__ import (absolute_import, print_function, division,
                        unicode_literals)


class GWFError(Exception):
    pass


class TargetExistsError(GWFError):

    def __init__(self, target):
        self.target = target

        message = (
            'Target "{}" already exists in workflow.'
        ).format(target.name)

        super(TargetExistsError, self).__init__(message)


class FileProvidedByMultipleTargetsError(GWFError):

    def __init__(self, path, target, other_target):
        self.path = path
        self.target = target
        self.other_target = other_target

        message = (
            'File "{}" provided by targets "{}" and "{}".'
        ).format(self.path, self.target, self.other_target)

        super(FileProvidedByMultipleTargetsError, self).__init__(message)


class FileRequiredButNotProvidedError(GWFError):

    def __init__(self, path, target):
        self.path = path
        self.target = target

        message = (
            'File "{}" is required by "{}", but does not exist and is not '
            'provided by any target in the workflow.'
        ).format(self.path, self.target)

        super(FileRequiredButNotProvidedError, self).__init__(message)


class CircularDependencyError(GWFError):

    def __init__(self, target):
        self.target = target

        message = (
            'Target {} depends on itself.'
        ).format(self.target)

        super(CircularDependencyError, self).__init__(message)


class WorkflowNotPreparedError(GWFError):
    pass


class IncludeWorkflowError(GWFError):
    pass
