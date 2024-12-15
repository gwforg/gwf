import time
from pathlib import Path

from gwf.cli import main


def test_run_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["run"])
    assert "Submitted target Target1" in result.output
    assert "Submitted target Target2" in result.output
    assert "Submitted target Target3" in result.output

    for _ in range(30):
        result = cli_runner.invoke(main, ["status", "-s", "shouldrun"])
        if result.output == "":
            return
        time.sleep(0.1)
    assert False


def test_run_dry_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["run", "--dry-run"])
    assert "Would submit Target1" in result.output
    assert "Would submit Target2" in result.output
    assert "Would submit Target3" in result.output


def test_run_partially_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["run", "Target2"])
    assert "Submitted target Target1" in result.output
    assert "Submitted target Target2" in result.output
    assert "Submitted target Target3" not in result.output
