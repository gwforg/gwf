"""
Various helper functions
"""

import string
import os
import os.path


def escape_file_name(filename):
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in filename if c in valid_chars)


def escape_job_name(job_name):
    valid_chars = "_%s%s" % (string.ascii_letters, string.digits)
    return ''.join(c for c in job_name if c in valid_chars)


_remembered_files = {}  # cache to avoid too much stat'ing
_remembered_timestamps = {}  # Use this to avoid too many stats that slows down the script


def file_exists(filename):
    if filename not in _remembered_files:
        _remembered_files[filename] = os.path.exists(filename)
    return _remembered_files[filename]


def get_file_timestamp(filename):
    if filename not in _remembered_timestamps:
        _remembered_timestamps[filename] = os.path.getmtime(filename)
    return _remembered_timestamps[filename]


def make_absolute_path(working_dir, filename):
    if os.path.isabs(filename):
        abspath = filename
    else:
        abspath = os.path.join(working_dir, filename)
    return os.path.normpath(abspath)


def _list(x):
    """Wrap x as a singleton in a list if it isn't a list already."""
    if hasattr(x, '__iter__'):
        return list(x)
    else:
        return [x]
