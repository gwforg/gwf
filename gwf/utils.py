import functools
import imp
import logging
import os.path
import re
import time
from contextlib import ContextDecorator

from gwf.exceptions import GWFError

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


def _split_import_path(path, default_obj):
    comps = path.rsplit(':')
    if len(comps) == 2:
        path, obj = comps
    elif len(comps) == 1:
        path, obj = comps[0], default_obj
    else:
        raise ValueError('Invalid path: "{}".'.format(path))

    basedir, filename = os.path.split(path)
    filename, ext = os.path.splitext(filename)
    return filename, basedir, obj


def import_object(path, default_obj='gwf'):
    if not os.path.isabs(path):
        path = os.path.abspath(os.path.join(os.getcwd(), path))

    filename, basedir, objname = _split_import_path(path, default_obj)
    fullpath = os.path.join(basedir, filename + '.py')

    if not os.path.exists(fullpath):
        raise GWFError(
            'The file "{}" does not exist.'.format(fullpath)
        )

    mod_loc = imp.find_module(filename, [basedir])
    mod = imp.load_module(filename, *mod_loc)

    try:
        return getattr(mod, objname)
    except AttributeError as e:
        logger.debug(e)
        raise GWFError(
            'Module "{}" does not declare the attribute "{}".'.format(
                filename,
                objname
            )
        )


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


def get_gwf_version():
    import pkg_resources
    return pkg_resources.get_distribution('gwf').version


def is_valid_name(candidate):
    return re.match(r'^[a-zA-Z_][a-zA-Z0-9._]*$', candidate) is not None


def merge(*args):
    first, *rest = args
    res = first.copy()
    for dct in rest:
        res.update(dct)
    return res


def safe_mkdir(path):
    try:
        os.makedirs(path)
    except OSError:
        pass
