import copy
import fnmatch
import functools
import importlib
import logging
import os.path
import re
import sys
import time
from contextlib import ContextDecorator

import click

from gwf.exceptions import GWFError, TargetDoesNotExistError

logger = logging.getLogger(__name__)


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
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.duration = self.end - self.start
        self.logger.debug(self.msg, self.duration * 1000)


def iter_inputs(targets):
    for target in targets:
        for path in target.inputs:
            yield target, path


def iter_outputs(targets):
    for target in targets:
        for path in target.outputs:
            yield target, path


def load_workflow(basedir, filename, objname):
    if not basedir:
        basedir = os.getcwd()
    fullpath = os.path.join(basedir, filename + '.py')

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
        raise GWFError('Module "{}" does not declare attribute "{}".'.format(filename, objname))


def get_file_timestamp(filename):
    """Return the modification time of `filename`.

    :param str filename: Path to a file.
    :return: Modification time of the file or `None` if the file does not exist.
    """
    try:
        return os.path.getmtime(filename)
    except OSError:
        return None


def dfs(root, dependencies):
    """Return the depth-first traversal path through a graph from `root`."""
    visited = set()
    path = []

    def dfs_inner(node):
        if node in visited:
            return

        visited.add(node)
        for dep in dependencies[node]:
            dfs_inner(dep)

        path.append(node)

    dfs_inner(root)
    return path


def is_valid_name(candidate):
    return re.match(r'^[a-zA-Z_][a-zA-Z0-9._]*$', candidate) is not None


def ensure_dir(path):
    """Create directory unless it already exists."""
    os.makedirs(path, exist_ok=True)


def parse_path(path, default_obj='gwf'):
    comps = path.rsplit(':')
    if len(comps) == 2:
        path, obj = comps
    elif len(comps) == 1:
        path, obj = comps[0], default_obj
    else:
        raise ValueError('Invalid path: "{}".'.format(path))

    basedir, filename = os.path.split(path)
    filename, _ = os.path.splitext(filename)
    return basedir, filename, obj


def match_targets(names, targets):
    matched_targets = set()
    for name in names:
        if '*' in name:
            for match in fnmatch.filter(targets.keys(), name):
                matched_targets.add(targets[match])
        elif name not in targets:
            raise TargetDoesNotExistError(name)
        else:
            matched_targets.add(targets[name])
    return matched_targets


class ColorFormatter(logging.Formatter):
    STYLING = {
        'WARNING': dict(fg='yellow'),
        'INFO': dict(fg='blue'),
        'DEBUG': dict(fg='white'),
        'ERROR': dict(fg='red', bold=True),
        'CRITICAL': dict(fg='magenta', bold=True),
    }

    def format(self, record):
        level = record.levelname
        color_record = copy.copy(record)
        if record.levelname in self.STYLING:
            styling = self.STYLING[level]
            color_record.levelname = click.style(record.levelname, **styling)
            color_record.msg = click.style(record.msg, **styling)
        return super().format(color_record)
