try:
    from os import fspath
except ImportError:
    def fspath(path):
        """For compatability with Python 3.5.

        Can be used instead of os.fspath() which is available in Python 3.6+.
        """
        if isinstance(path, (str, bytes)):
            return path
        path = path.__fspath__()
        if isinstance(path, (str, bytes)):
            return path
        raise TypeError('not a path')
