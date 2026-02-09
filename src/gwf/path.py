import os
from os import PathLike
from pathlib import Path


def _standardize_path(path: str | PathLike | Path) -> str:
    if isinstance(path, Path):
        return str(path)
    return os.fspath(path)


class BasePath:
    def __init__(self, path: str | PathLike | Path):
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
    def __init__(self, path: str | PathLike | Path):
        super().__init__(path)


class TemporaryPath(BasePath):
    def __init__(self, path: str | PathLike | Path):
        super().__init__(path)


def temp(path: str | PathLike | Path) -> TemporaryPath:
    """Mark a path as temporary.

    Temporary files will be automatically cleaned up after the workflow completes.

    Args:
        path: The path to mark as temporary. Can be a string, PathLike, or Path object.

    Returns:
        A TemporaryPath instance wrapping the given path.

    Example:
        >>> outputs = {"intermediate": temp("data/intermediate.txt")}
    """
    return TemporaryPath(path)


def protect(path: str | PathLike | Path) -> ProtectedPath:
    """Mark a path as protected from cleanup.

    Protected files will not be removed during cleanup operations, even if they
    are outputs of non-endpoint targets.

    Args:
        path: The path to mark as protected. Can be a string, PathLike, or Path object.

    Returns:
        A ProtectedPath instance wrapping the given path.

    Example:
        >>> outputs = {"important": protect("results/important_data.txt")}
    """
    return ProtectedPath(path)
