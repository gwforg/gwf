import copy
import functools
import importlib
import itertools
import json
import logging
import os
import os.path
import re
import socket
import sys
import time
from collections import UserDict
from contextlib import ContextDecorator
from functools import wraps
from urllib.request import urlopen

import click

from gwf.exceptions import GWFError

UPDATE_CHECK_URL = "https://pypi.org/pypi/gwf/json"
UPDATE_CHECK_FILE = ".gwf/update"
UPDATE_CHECK_WAIT = 24 * 60 * 60


logger = logging.getLogger(__name__)


def is_valid_name(candidate):
    """Check whether `candidate` is a valid name for a target or workflow."""
    return re.match(r"^[a-zA-Z_][a-zA-Z0-9._]*$", candidate) is not None


def cache(obj):
    _cache = obj._cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        if args not in _cache:
            _cache[args] = obj(*args, **kwargs)
        return _cache[args]

    return memoizer


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


def parse_path(path, default_obj="gwf", default_file="workflow.py"):
    comps = path.rsplit(":")
    if len(comps) == 2:
        path, obj = comps[0] or default_file, comps[1] or default_obj
    elif len(comps) == 1:
        path, obj = comps[0], default_obj
    else:
        raise ValueError('Invalid path: "{}".'.format(path))

    basedir, filename = os.path.split(path)
    filename, _ = os.path.splitext(filename)
    return basedir, filename, obj


def load_workflow(basedir, filename, objname):
    if not basedir:
        basedir = os.getcwd()
    fullpath = os.path.join(basedir, filename + ".py")

    if not os.path.exists(fullpath):
        raise GWFError('The file "{}" does not exist.'.format(fullpath))

    sys.path.insert(0, os.path.join(os.getcwd(), basedir))
    spec = importlib.util.spec_from_file_location(filename, fullpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sys.path.pop(0)

    try:
        return getattr(mod, objname)
    except AttributeError:
        raise GWFError(
            'Module "{}" does not declare attribute "{}".'.format(filename, objname)
        )


class PersistableDict(UserDict):
    """A dictionary which can persist itself to JSON."""

    def __init__(self, path):
        super().__init__()

        self.path = path
        try:
            with open(self.path) as fileobj:
                self.data.update(json.load(fileobj))
        except (OSError, ValueError):
            # Catch ValueError for compatibility with Python 3.4.2. I haven't been
            # able to figure out what is different between 3.4.2 and 3.5 that
            # causes this. Essentially, 3.4.2 raises a ValueError saying that it
            # cannot parse the empty string instead of raising an OSError
            # (FileNotFoundError does not exist in 3.4.2) saying that the file does
            # not exist.
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


def get_latest_version():
    """Return the latest version available.

    Will return a string containing the version of the latest release.

    Contacting the server to check for updates will happen at most once per
    day. If the last check happened within 24 hours, `None` will be returned.

    :return: A string containing the version number or `None`.
    """
    try:
        last_check = os.stat(UPDATE_CHECK_FILE).st_mtime
    except FileNotFoundError:
        last_check = 0

    current_time = time.time()
    if current_time < last_check + UPDATE_CHECK_WAIT:
        logger.debug("Skipping check for updates.")
        return None

    logger.debug("Checking for updates.")
    touchfile(UPDATE_CHECK_FILE)
    try:
        with urlopen(UPDATE_CHECK_URL, timeout=1) as resp:
            data = json.load(resp)
            latest_version = data["info"]["version"]
            return latest_version
    except socket.timeout:
        logger.debug("Connect to version server timed out.")
        return None


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

                    delay = (2 ** retries) // 2
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
