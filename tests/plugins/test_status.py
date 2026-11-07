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


def test_status_summary_when_no_targets_match(
    cli_runner, local_backend, simple_workflow
):
    result = cli_runner.invoke(main, ["status", "--format", "summary", "NoSuchTarget"])

    assert result.exit_code == 0, result.output
    assert "completed" in result.output


def test_status_aligns_status_column_to_longest_target_name(
    cli_runner, local_backend, tmpdir
):
    tmpdir.join("workflow.py").write(
        """from gwf import Workflow

gwf = Workflow()
gwf.target('A', inputs=[], outputs=['a.txt']) << 'touch a.txt'
gwf.target('AMuchLongerTargetName', inputs=[], outputs=['b.txt']) << 'touch b.txt'
"""
    )

    with tmpdir.as_cwd():
        result = cli_runner.invoke(main, ["status"])

    assert result.exit_code == 0, result.output
    lines = [line for line in result.output.splitlines() if "shouldrun" in line]
    assert len(lines) == 2
    assert len({line.index("shouldrun") for line in lines}) == 1
