import os
from os import PathLike
from pathlib import Path

TEMP_PATH_REGISTRY = set()


def _standardize_path(path: str | PathLike | Path) -> str:
    if isinstance(path, Path):
        return str(path)
    return os.fspath(path)


def _clear_temp_registry():
    TEMP_PATH_REGISTRY.clear()


def temp(path: str | PathLike | Path) -> str | PathLike | Path:
    TEMP_PATH_REGISTRY.add(_standardize_path(path))
    return path


def is_temp(path) -> bool:
    return _standardize_path(path) in TEMP_PATH_REGISTRY
