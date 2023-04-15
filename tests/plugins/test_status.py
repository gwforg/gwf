from pathlib import Path

from gwf.cli import main


def test_status_shows_all_targets(cli_runner, local_backend, simple_workflow):
    result = cli_runner.invoke(main, ["status"])
    assert "Target1" in result.output
    assert "Target2" in result.output
    assert "Target3" in result.output


def test_status_shows_one_named_target(cli_runner, local_backend, simple_workflow):
    result = cli_runner.invoke(main, ["status", "Target1"])
    assert "Target1" in result.output
    assert "Target2" not in result.output
    assert "Target3" not in result.output


def test_status_shows_two_named_targets(cli_runner, local_backend, simple_workflow):
    result = cli_runner.invoke(main, ["status", "Target1", "Target2"])
    assert "Target2" in result.output
    assert "Target1" in result.output
    assert "Target3" not in result.output


def test_status_shows_only_endpoint_targets(cli_runner, local_backend, simple_workflow):
    result = cli_runner.invoke(main, ["status", "--endpoints"])
    assert "Target1" not in result.output
    assert "Target2" in result.output
    assert "Target3" in result.output


def test_status_shows_only_completed_targets(
    cli_runner, local_backend, linear_workflow
):
    result = cli_runner.invoke(main, ["config", "set", "use_spec_hashes", "no"])

    Path("a.txt").touch()
    Path("b.txt").touch()
    Path("c.txt").touch()

    result = cli_runner.invoke(main, ["status", "-s", "completed"])
    assert "Target1" in result.output
    assert "Target2" in result.output
    assert "Target3" not in result.output
