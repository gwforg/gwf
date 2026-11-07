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


GROUPED_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('A', inputs=[], outputs=['a.txt'], group='first') << 'touch a.txt'
gwf.target('B', inputs=['a.txt'], outputs=['b.txt'], group='second') << 'touch b.txt'
"""


def test_run_group_only_submits_targets_in_group(cli_runner, local_backend, tmpdir):
    tmpdir.join("workflow.py").write(GROUPED_WORKFLOW)

    with tmpdir.as_cwd():
        result = cli_runner.invoke(main, ["run", "--dry-run", "-g", "first"])

    assert result.exit_code == 0, result.output
    assert "Would submit A" in result.output
    assert "Would submit B" not in result.output


def test_run_status_only_submits_targets_with_status(cli_runner, local_backend, tmpdir):
    tmpdir.join("workflow.py").write(GROUPED_WORKFLOW)

    with tmpdir.as_cwd():
        Path("a.txt").touch()
        result = cli_runner.invoke(main, ["run", "--dry-run", "-s", "completed"])

    assert result.exit_code == 0, result.output
    assert "Would submit" not in result.output


def test_run_names_and_group_must_both_match(cli_runner, local_backend, tmpdir):
    tmpdir.join("workflow.py").write(GROUPED_WORKFLOW)

    with tmpdir.as_cwd():
        result = cli_runner.invoke(main, ["run", "--dry-run", "-g", "first", "B"])

    assert result.exit_code == 0, result.output
    assert "Would submit" not in result.output
