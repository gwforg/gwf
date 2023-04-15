import copy
import itertools
import json
import logging
import os
import os.path
import re
import sys
import time
from collections import UserDict
from contextlib import ContextDecorator
from functools import wraps
from pathlib import Path

if sys.version_info < (3, 10):
    from importlib_metadata import entry_points  # noqa: E401
else:
    from importlib.metadata import entry_points  # noqa: F401

import click

from gwf.exceptions import GWFError

logger = logging.getLogger(__name__)


def is_valid_name(candidate):
    """Check whether `candidate` is a valid name for a target or workflow."""
    return re.match(r"^[a-zA-Z_][a-zA-Z0-9._]*$", candidate) is not None


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


class PersistableDict(UserDict):
    """A dictionary which can persist itself to JSON."""

    def __init__(self, path):
        super().__init__()

        self.path = path
        try:
            with open(self.path) as fileobj:
                self.data.update(json.load(fileobj))
        except OSError:
            pass

    def persist(self):
        with open(self.path + ".new", "w") as fileobj:
            json.dump(self.data, fileobj)
            fileobj.flush()
            os.fsync(fileobj.fileno())
            fileobj.close()
        os.rename(self.path + ".new", self.path)


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


def touchfile(path):
    with open(path, "a"):
        os.utime(path, None)


class RetryError(Exception):
    """Raised when max_retries has been exceeded."""


def retry(on_exc, max_retries=3, callback=None):
    """Retry a function with exponentially increasing delay.

    This will retry the decorated function up to `max_retries` times. A retry
    will only be attempted if the exception raised by the wrapped function is
    in `on_exc`.

    If `callback(retries, delay)` is given, it must be a callable the number of
    retries so far as the first argument and the current delay as the second
    argument. Its return value is ignored.
    """

    def retry_decorator(func):
        @wraps(func)
        def func_wrapper(*args, **kwargs):
            last_exc = None
            for retries in itertools.count():  # pragma: no cover
                if retries >= max_retries:
                    raise RetryError(func.__name__) from last_exc

                try:
                    return func(*args, **kwargs)
                except on_exc as exc:
                    last_exc = exc

                    delay = (2**retries) // 2
                    if callback is not None:
                        callback(retries, delay)

                    logger.exception(exc)
                    logger.warning(
                        "Call to %s failed, retrying in %d seconds.",
                        func.__name__,
                        delay,
                    )
                    time.sleep(delay)

        return func_wrapper

    return retry_decorator


retry.RetryError = RetryError


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
                raise GWFError(f"The file {path} could not be found")
            current_dir = current_dir.parent
            workflow_path = current_dir.joinpath(path)
    return current_dir, workflow_path, obj
