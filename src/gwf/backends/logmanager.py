import io
import os
import os.path
import threading

from ..utils import ensure_dir, redirect_exception
from .exceptions import LogError


class OutputStream(threading.Thread):
    def __init__(self):
        super().__init__()
        self.done = False
        self.buffer = io.StringIO()
        self.read_fd, self.write_fd = os.pipe()
        self.reader = os.fdopen(self.read_fd)
        self.start()

    def fileno(self):
        return self.write_fd

    def run(self):
        while not self.done:
            self.buffer.write(self.reader.read())
        self.reader.close()

    def close(self):
        self.done = True
        os.close(self.write_fd)


class LogManager:
    """Base class for log managers.

    A log manager must as a minimum implement the :func:`open_stdout(target)`
    and :func:`open_stderr(target)` methods.

    Log managers are used by backends to specify how log files are stored and
    retrieved. Most backends will want to use the :class:`FileLogManager` which
    stores log files in the `.gwf/logs` directory. For testing of backends, the
    :class:`MemoryLogManager` can be useful.
    """

    def open_stdout(self, target, mode="r"):
        raise NotImplementedError()

    def open_stderr(self, target, mode="r"):
        raise NotImplementedError()

    def remove_stdout(self, target):
        raise NotImplementedError()

    def remove_stderr(self, target):
        raise NotImplementedError()

    def list(self):
        raise NotImplementedError()


class MemoryLogManager(LogManager):
    """A memory-based log manager.

    This log manager stores logs in memory.
    """

    def __init__(self):
        self._stdout_logs = {}
        self._stderr_logs = {}

    @redirect_exception(KeyError, LogError)
    def open_stdout(self, target, mode="r"):
        if target not in self._stdout_logs and mode == "w":
            self._stdout_logs[target] = OutputStream()
        return self._stdout_logs[target]

    @redirect_exception(KeyError, LogError)
    def open_stderr(self, target, mode="r"):
        if target not in self._stderr_logs and mode == "w":
            self._stderr_logs[target] = OutputStream()
        return self._stderr_logs[target]

    @redirect_exception(KeyError, LogError)
    def remove_stdout(self, target):
        del self._stdout_logs[target]

    @redirect_exception(KeyError, LogError)
    def remove_stderr(self, target):
        del self._stderr_logs[target]

    def list(self):
        return set(self._stderr_logs.keys()).union(self._stdout_logs.keys())


class FileLogManager(LogManager):
    """A file-based log manager.

    This log manager stores logs on disk in the `log_dir` directory (which
    defaults to `.gwf/logs`).
    """

    log_dir = ".gwf/logs"

    def __init__(self):
        super().__init__()
        ensure_dir(FileLogManager.log_dir)

    @staticmethod
    def _get_log_path(target_name, extension):
        """Return path for log file for a given target name.

        If `extension` is `stdout`, then the path of the file containing
        standard output of the target will be returned. If `stderr`, the path of
        the file containing standard error of the target will be returned.

        :arg target gwf.Target:
            Target to return log path for.
        :arg extension str:
            Must be either `stdout` or `stderr`.
        """
        log_file = "{}.{}".format(target_name, extension)
        return os.path.join(FileLogManager.log_dir, log_file)

    def stdout_path(self, target):
        """Return path of the log file containing standard output for target."""
        if hasattr(target, "name"):
            target = target.name
        return self._get_log_path(target, "stdout")

    def stderr_path(self, target):
        """Return path of the log file containing standard error for target."""
        if hasattr(target, "name"):
            target = target.name
        return self._get_log_path(target, "stderr")

    @redirect_exception(FileNotFoundError, LogError)
    def open_stdout(self, target, mode="r"):
        """Return file handle to the standard output log file for target.

        :raises LogError: If the log could not be found.
        """
        if hasattr(target, "name"):
            target = target.name
        return open(self.stdout_path(target), mode)

    @redirect_exception(FileNotFoundError, LogError)
    def open_stderr(self, target, mode="r"):
        """Return file handle to standard error log file for target.

        :raises LogError: If the log could not be found.
        """
        if hasattr(target, "name"):
            target = target.name
        return open(self.stderr_path(target), mode)

    @redirect_exception(OSError, LogError)
    def remove_stdout(self, target):
        if hasattr(target, "name"):
            target = target.name
        os.remove(self.stdout_path(target))

    @redirect_exception(OSError, LogError)
    def remove_stderr(self, target):
        if hasattr(target, "name"):
            target = target.name
        os.remove(self.stderr_path(target))

    def list(self):
        return (
            os.path.splitext(log_name)[0]
            for log_name in os.listdir(FileLogManager.log_dir)
        )
