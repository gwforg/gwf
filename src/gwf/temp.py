import os
from os import PathLike
from pathlib import Path


def _standardize_path(path: str | PathLike | Path) -> str:
    if isinstance(path, Path):
        return str(path)
    return os.fspath(path)


class TempRegistry(set):
    """A registry to keep track of temporary paths."""

    def __contains__(self, o):
        return super().__contains__(_standardize_path(o))

    def register(self, path) -> None:
        """Register a temporary path."""
        self.add(_standardize_path(path))

    def unregister(self, path) -> None:
        """Unregister a temporary path."""
        self.discard(_standardize_path(path))

    def update(self, paths: dict) -> None:
        """Update registered paths based on a mapping from current to new paths."""
        for current_path, new_path in paths.items():
            self.discard(_standardize_path(current_path))
            self.add(_standardize_path(new_path))


_temp_registry = TempRegistry()


def temp(path: str | PathLike | Path) -> str:
    """Register a path as temporary."""
    _temp_registry.register(path)
    return path
