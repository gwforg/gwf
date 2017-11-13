import functools
import importlib
import json
import logging
import os.path
import sys
import time
from collections import UserDict
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


class LazyDict(dict):
    """A dict which lazily computes values for keys using `valfunc`.

    When accessing an key in the dict, it will check whether the key exists. If it does, the value is returned
    immediately. If not, `valfunc` will be called on the key and the return value will be assigned as the value of the
    key. For example::

        >>> d = LazyDict(valfunc=lambda k: k + 1)
        >>> 0 in d
        False
        >>> d[0]
        1
        >>> d[100]
        101
    """
    def __init__(self, valfunc, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.valfunc = valfunc

    def __getitem__(self, item):
        if not super().__contains__(item):
            super().__setitem__(item, self.valfunc(item))
        return super().__getitem__(item)


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
        with open(self.path + '.new', 'w') as fileobj:
            json.dump(self.data, fileobj)
            fileobj.flush()
            os.fsync(fileobj.fileno())
            fileobj.close()
        os.rename(self.path + '.new', self.path)
