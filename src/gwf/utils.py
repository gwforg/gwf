import copy
import importlib
import logging
import os
import os.path
import re
import sys
import time
from contextlib import ContextDecorator
from functools import wraps
from pathlib import Path

if sys.version_info < (3, 8):
    from importlib_metadata import entry_points as _entry_points  # noqa: E401
else:
    from importlib.metadata import entry_points as _entry_points  # noqa: F401

import click

from gwf.exceptions import GWFError

logger = logging.getLogger(__name__)


def entry_points(group):
    if sys.version_info < (3, 12):
        return _entry_points()[group]
    else:
        return _entry_points(group=group)


def is_valid_name(candidate):
    """Check whether `candidate` is a valid name for a target or workflow."""
    return re.match(r"^[a-zA-Z_][a-zA-Z0-9._]*$", candidate) is not None


def chain(*dcts):
    new = {}
    for dct in dcts:
        new.update(dct)
    return new


class timer(ContextDecorator):
    def __init__(self, msg, logger=None):
        self.msg = msg
        self.logger = logger or logging.getLogger(__name__)

    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *args):
        self.end = time.perf_counter()
        self.duration = self.end - self.start
        self.logger.debug(self.msg, self.duration)


def ensure_dir(path):
    """Create directory unless it already exists."""
    os.makedirs(path, exist_ok=True)


class ColorFormatter(logging.Formatter):
    STYLING = {
        "WARNING": dict(fg="yellow"),
        "INFO": dict(fg="blue"),
        "DEBUG": dict(fg="cyan"),
        "ERROR": dict(fg="red", bold=True),
        "CRITICAL": dict(fg="magenta", bold=True),
    }

    def format(self, record):
        level = record.levelname
        color_record = copy.copy(record)
        if record.levelname in self.STYLING:
            styling = self.STYLING[level]
            padded_level_name = "{:<10}".format(record.levelname.lower())
            color_record.levelname = click.style(padded_level_name, **styling)
            color_record.name = click.style(record.name, **styling)
            color_record.msg = record.msg
        return super().format(color_record)


def redirect_exception(old_exc, new_exc):
    """Redirect one exception type to another."""

    def wrapper(func):
        @wraps(func)
        def inner_wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except old_exc as e:
                raise new_exc from e

        return inner_wrapper

    return wrapper


def ensure_trailing_newline(s):
    if not s:
        return "\n"
    return s if s[-1] == "\n" else s + "\n"


def find_workflow(path_spec):
    path, _, obj = path_spec.partition(":")
    obj = obj or "gwf"

    path = Path(path)
    current_dir = Path.cwd()
    workflow_path = current_dir.joinpath(path)
    if not path.is_absolute():
        while True:
            logger.debug("Looking for workflow file %s in %s", path, current_dir)
            if workflow_path.exists():
                break
            if current_dir == Path(current_dir.anchor):
                raise FileNotFoundError(f"The file {path} could not be found")
            current_dir = current_dir.parent
            workflow_path = current_dir.joinpath(path)
    return workflow_path, obj


def load_workflow(path: Path, obj: str):
    logger.debug("Loading workflow from %s:%s", path, obj)
    module_name, _ = os.path.splitext(path.name)
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(module_name, path)
    assert spec is not None, "Could not load workflow file"
    module = importlib.util.module_from_spec(spec)
    assert module is not None, "Could not load module from workflow file"
    spec.loader.exec_module(module)
    if not hasattr(module, obj):
        raise GWFError(f"The module '{path.name}' does not have attribute '{obj}'")
    return getattr(module, obj)
