from pathlib import Path

from gwf.cli import main


def test_logs_shows_log_file_from_target(cli_runner, linear_workflow):
    Path.cwd().joinpath(".gwf", "logs").mkdir(parents=True)
    Path.cwd().joinpath(".gwf", "logs", "Target1.stdout").write_text("hello")
    Path.cwd().joinpath(".gwf", "logs", "Target1.stderr").write_text("world")

    result = cli_runner.invoke(main, ["logs", "Target1"])
    assert "hello" in result.output

    result = cli_runner.invoke(main, ["logs", "-e", "Target1"])
    assert "world" in result.output
