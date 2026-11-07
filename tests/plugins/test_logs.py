from pathlib import Path

from gwf.cli import main
from gwf.log_storage import (
    get_log_paths,
    get_logs_dir,
    prepare_log_storage_for_target,
)


def test_logs_shows_log_file_from_target(cli_runner, linear_workflow):
    logs_dir = get_logs_dir(Path.cwd())
    logs_dir.mkdir(parents=True)
    logs_dir = logs_dir.joinpath("a", "2")
    logs_dir.mkdir(parents=True)
    logs_dir.joinpath("Target1.stdout").write_text("hello")
    logs_dir.joinpath("Target1.stderr").write_text("world")

    result = cli_runner.invoke(main, ["logs", "Target1"])
    assert "hello" in result.output

    result = cli_runner.invoke(main, ["logs", "-e", "Target1"])
    assert "world" in result.output


def test_logs_migrated_to_new_hierarchy(cli_runner, linear_workflow):
    logs_dir = get_logs_dir(Path.cwd())
    logs_dir.mkdir(parents=True)

    # Write logs to in old location.
    logs_dir.joinpath("Target1.stdout").write_text("hello")
    logs_dir.joinpath("Target1.stderr").write_text("world")

    cli_runner.invoke(main, ["logs", "Target1"])

    # The logs should now exist in the new, hash-based location.    print(logs_dir)
    assert logs_dir.joinpath("a", "2", "Target1.stdout").exists()
    assert logs_dir.joinpath("a", "2", "Target1.stderr").exists()


def _write_log(name, stdout=b"", stderr=b""):
    prepare_log_storage_for_target(Path.cwd(), name)
    stdout_path, stderr_path = get_log_paths(Path.cwd(), name)
    stdout_path.write_bytes(stdout)
    stderr_path.write_bytes(stderr)


def test_logs_without_pager_shows_log_contents(cli_runner, linear_workflow):
    _write_log("Target1", stdout=b"hello\nworld\n")

    result = cli_runner.invoke(main, ["logs", "--no-pager", "Target1"])

    assert result.exit_code == 0, result.output
    assert result.output == "hello\nworld\n"


def test_logs_reports_missing_log(cli_runner, linear_workflow):
    result = cli_runner.invoke(main, ["logs", "--no-pager", "-e", "Target1"])

    assert result.exit_code == 1
    assert "No standard error log found for target Target1" in result.output
    assert not isinstance(result.exception, FileNotFoundError)


def test_logs_tolerates_invalid_utf8(cli_runner, linear_workflow):
    _write_log("Target1", stdout=b"before \xff after\n")

    result = cli_runner.invoke(main, ["logs", "--no-pager", "Target1"])

    assert result.exit_code == 0, result.output
    assert "before" in result.output
    assert "after" in result.output
