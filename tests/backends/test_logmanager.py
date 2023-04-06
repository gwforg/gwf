import pytest

from gwf import Target
from gwf.backends.exceptions import LogError
from gwf.backends.logmanager import FileLogManager


@pytest.fixture
def log_manager(tmpdir):
    with tmpdir.as_cwd():
        yield FileLogManager()


@pytest.fixture
def empty_target():
    return Target(
        name="TestTarget", inputs=[], outputs=[], options={}, working_dir="/some/dir"
    )


def test_file_log_manager_creates_logs_dir(log_manager, tmpdir):
    assert tmpdir.join(".gwf", "logs").exists()


def test_file_log_manager_path(log_manager, empty_target):
    assert log_manager.stderr_path(empty_target) == ".gwf/logs/TestTarget.stderr"
    assert log_manager.stdout_path(empty_target) == ".gwf/logs/TestTarget.stdout"


def test_file_log_manager_open(log_manager, tmpdir, empty_target):
    with pytest.raises(LogError):
        log_manager.open_stdout(empty_target)

    with pytest.raises(LogError):
        log_manager.open_stderr(empty_target)

    stdout_file = tmpdir.join(".gwf", "logs", "TestTarget.stdout")
    stdout_file.write("This is standard output...")

    stderr_file = tmpdir.join(".gwf", "logs", "TestTarget.stderr")
    stderr_file.write("This is standard error...")

    assert log_manager.open_stdout(empty_target).read() == "This is standard output..."
    assert log_manager.open_stderr(empty_target).read() == "This is standard error..."


def test_file_log_manager_remove(log_manager, tmpdir, empty_target):
    stdout_file = tmpdir.join(".gwf", "logs", "TestTarget.stdout").ensure()
    stderr_file = tmpdir.join(".gwf", "logs", "TestTarget.stderr").ensure()

    log_manager.remove_stdout(empty_target)
    assert not stdout_file.exists()

    log_manager.remove_stderr(empty_target)
    assert not stderr_file.exists()

    with pytest.raises(LogError):
        log_manager.remove_stdout(empty_target)

    with pytest.raises(LogError):
        log_manager.remove_stderr(empty_target)
