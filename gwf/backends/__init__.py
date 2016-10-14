import os.path

from ..exceptions import NoLogFoundError
from .base import Backend


class FileLogsMixin:

    def logs(self, target, stderr=False):
        try:
            if stderr:
                return self.open_stderr(target)
            return self.open_stdout(target)
        except OSError as e:
            raise NoLogFoundError('Could not find logs.') from e

    @staticmethod
    def log_dir(target):
        """Path to directory containing logs for `target`."""
        return os.path.join(target.working_dir, '.gwf', 'logs')

    @staticmethod
    def open_stdout(target, mode='r'):
        path = os.path.join(
            FileLogsMixin.log_dir(target),
            '{}.stdout'.format(target.name)
        )
        return open(path, mode)

    @staticmethod
    def open_stderr(target, mode='r'):
        path = os.path.join(
            FileLogsMixin.log_dir(target),
            '{}.stderr'.format(target.name)
        )
        return open(path, mode)

__all__ = ('Backend', 'FileLogsMixin',)
