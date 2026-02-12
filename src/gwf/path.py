import os
from os import PathLike
from pathlib import Path


def _standardize_path(path) -> str:
    if isinstance(path, Path):
        return str(path)
    return os.fspath(path)


class BasePath:
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return _standardize_path(self.path)

    def __fspath__(self):
        return _standardize_path(self.path)

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))


class ProtectedPath(BasePath):
    def __init__(self, path):
        super().__init__(path)


class TemporaryPath(BasePath):
    def __init__(self, path):
        super().__init__(path)


def temp(path) -> TemporaryPath:
    """Mark a path as temporary.

    Temporary files are expected to be intermediate files that are not needed after the workflow
    has completed. Temporary files can be removed after the workflow finishes without affecting
    the scheduling states.

    :param path:
        The path to mark as temporary. Can be a string, PathLike, or Path object.

    :returns:
        A TemporaryPath instance wrapping the given path.

    Example:
        >>> outputs = {"intermediate": temp("data/intermediate.txt")}
    """
    return TemporaryPath(path)


def protect(path) -> ProtectedPath:
    """Mark a path as protected from cleanup.

    Protected files will not be removed during cleanup operations, even if they are outputs of
    non-endpoint targets.

    :param path:
        The path to mark as protected. Can be a string, PathLike, or Path object.

    :returns:
        A ProtectedPath instance wrapping the given path.

    Example:
        >>> outputs = {"important": protect("results/important_data.txt")}
    """
    return ProtectedPath(path)
