import pytest

from gwf import Target
from gwf.backends.exceptions import LogNotFoundError
from gwf.backends.logmanager import FileLogManager


def test_file_log_manager(tmpdir):
    target = Target.empty('TestTarget')

    log_manager = FileLogManager()
    assert log_manager.stderr_path(target) == '.gwf/logs/TestTarget.stderr'
    assert log_manager.stdout_path(target) == '.gwf/logs/TestTarget.stdout'

    log_dir = tmpdir.join('.gwf', 'logs')
    log_dir.ensure_dir()

    with tmpdir.as_cwd():
        with pytest.raises(LogNotFoundError):
            log_manager.open_stdout(target)

        with pytest.raises(LogNotFoundError):
            log_manager.open_stderr(target)

        stdout_file = log_dir.join('TestTarget.stdout')
        stdout_file.write('This is standard output...')

        stderr_file = log_dir.join('TestTarget.stderr')
        stderr_file.write('This is standard error...')

        assert log_manager.open_stdout(target).read() == 'This is standard output...'
        assert log_manager.open_stderr(target).read() == 'This is standard error...'
