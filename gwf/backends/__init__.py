import os.path

from ..exceptions import NoLogFoundError
from ..utils import safe_mkdir
from .base import Backend


class FileLogsMixin:

    def logs(self, target, stderr=False):
        try:
            if stderr:
                return self.open_stderr(target)
            return self.open_stdout(target)
        except OSError as e:
            raise NoLogFoundError() from e

    @staticmethod
    def log_dir(target):
        """Path to directory containing logs for `target`."""
        return os.path.join(target.working_dir, '.gwf', 'logs')

    @staticmethod
    def _log_path(target, extension):
        log_dir = FileLogsMixin.log_dir(target)
        safe_mkdir(log_dir)
        return os.path.join(log_dir, '{}.{}'.format(target.name, extension))

    @staticmethod
    def stdout_path(target):
        return FileLogsMixin._log_path(target, extension='stdout')

    def stderr_path(target):
        return FileLogsMixin._log_path(target, extension='stderr')

    @staticmethod
    def open_stdout(target, mode='r'):
        return open(FileLogsMixin.stdout_path(target), mode)

    @staticmethod
    def open_stderr(target, mode='r'):
        return open(FileLogsMixin.stderr_path(target), mode)

__all__ = ('Backend', 'FileLogsMixin',)
