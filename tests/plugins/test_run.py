from pathlib import Path

from gwf.cli import main


def test_run_submits_targets(cli_runner, local_backend, linear_workflow, wait_for_path):
    Path("a.txt").touch()

    result = cli_runner.invoke(main, ["run"])
    assert "Submitted target Target1" in result.output
    assert "Submitted target Target2" in result.output
    assert "Submitted target Target3" in result.output

    wait_for_path("d.txt")


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


def test_run_can_ignore_unneeded_dependency_outputs(cli_runner, local_backend, tmpdir):
    tmpdir.join("workflow.py").write(
        """from gwf import Workflow

gwf = Workflow()
gwf.target('Producer', inputs=[], outputs=['x.txt', 'y.txt']) << 'touch x.txt y.txt'
gwf.target('NeedsX', inputs=['x.txt'], outputs=['b.txt']) << 'touch b.txt'
gwf.target('NeedsY', inputs=['y.txt'], outputs=['c.txt']) << 'touch c.txt'
"""
    )

    with tmpdir.as_cwd():
        Path("x.txt").touch()
        result = cli_runner.invoke(
            main,
            ["config", "set", "require_all_outputs", "false"],
        )
        assert result.exit_code == 0, result.output

        result = cli_runner.invoke(main, ["run", "--dry-run", "NeedsX"])

    assert result.exit_code == 0, result.output
    assert "Would submit Producer" not in result.output
    assert "Would submit NeedsX" in result.output
