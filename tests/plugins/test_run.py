from pathlib import Path

from gwf.cli import main


def test_run_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["-b", "local", "run"])
    print(result.exception)
    assert "Submitting target Target1" in result.output
    assert "Submitting target Target2" in result.output
    assert "Submitting target Target3" in result.output


def test_run_dry_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["-b", "local", "run", "--dry-run"])
    assert "Would submit Target1" in result.output
    assert "Would submit Target2" in result.output
    assert "Would submit Target3" in result.output


def test_run_partially_submits_targets(cli_runner, local_backend, linear_workflow):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["-b", "local", "run", "Target2"])
    assert "Submitting target Target1" in result.output
    assert "Submitting target Target2" in result.output
    assert "Submitting target Target3" not in result.output
