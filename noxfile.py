import nox

import os
import shutil
from fnmatch import fnmatch


PURGE_PATTERNS = [
    "__pycache__/*",
    "docs/_build/*",
    ".coverage/*",
    ".egg-info/*",
    ".egg/*",
    ".eggs/*",
    ".gwf/*",
    ".gwfconf.json",
    ".pytest_cache",
    "./conda-bld",
    "dist",
    "*.pyc",
]


def matches_pattern(path):
    for pattern in PURGE_PATTERNS:
        if fnmatch(path, pattern):
            return True
    return False


@nox.session
def clean(session):
    def delete_recursively(path):
        if os.path.isdir(path):
            for child in os.scandir(path):
                delete_recursively(child.path)
            print("Deleting directory", path)
            os.rmdir(path)
        else:
            print("Deleting file", path)
            os.remove(path)

    for root, dirs, files in os.walk(".", topdown=False):
        for name in dirs + files:
            if matches_pattern(name):
                delete_recursively(os.path.join(root, name))


@nox.session
def build(session):
    session.install("flit")
    session.run("flit", "build", "--no-setup-py", "--format", "wheel")
