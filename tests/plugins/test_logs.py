from pathlib import Path

from gwf.cli import main
from gwf.log_storage import get_logs_dir


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

    result = cli_runner.invoke(main, ["logs", "Target1"])

    # The logs should now exist in the new, hash-based location.    print(logs_dir)
    assert logs_dir.joinpath("a", "2", "Target1.stdout").exists()
    assert logs_dir.joinpath("a", "2", "Target1.stderr").exists()
