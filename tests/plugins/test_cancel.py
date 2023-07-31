import pytest

from gwf.cli import main


@pytest.fixture
def long_running_workflow(tmpdir):
    workflow_file = tmpdir.join("workflow.py")
    workflow_file.write(
        """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << 'touch a.txt; sleep 3'
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt']) << 'sleep 3; touch b.txt'
"""
    )
    with tmpdir.as_cwd():
        yield tmpdir


def test_cancel_one_target(cli_runner, local_backend, long_running_workflow):
    result = cli_runner.invoke(main, ["run", "Target1"])
    result = cli_runner.invoke(main, ["cancel", "Target1"])
    assert result.output == "Cancelling target Target1\n"
    result = cli_runner.invoke(main, ["status", "Target1"])
    assert "cancelled" in result.output


def test_cancel_no_targets_specified_should_ask_for_confirmation_and_cancel_all_if_approved(
    cli_runner, local_backend, long_running_workflow
):
    result = cli_runner.invoke(main, ["cancel"], input="y")
    lines = result.output.strip().split("\n")
    assert "Cancelling target Target1" in lines
    assert "Cancelling target Target2" in lines


def test_cancel_no_targets_specified_should_ask_for_confirmation_and_abort_if_not_approved(
    cli_runner, local_backend, long_running_workflow
):
    result = cli_runner.invoke(main, ["cancel"], input="N")
    assert "Aborted!\n" in result.output
