import pytest

from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=[]) << "echo hello world"
gwf.target('Target2', inputs=[], outputs=[]) << "echo world hello"
"""


@pytest.fixture
def simple_workflow(tmpdir):
    workflow_file = tmpdir.join("workflow.py")
    workflow_file.write(SIMPLE_WORKFLOW)
    return tmpdir


@pytest.fixture(autouse=True)
def setup(simple_workflow):
    with simple_workflow.as_cwd():
        yield


def test_run_all_targets(cli_runner, mocker):
    mock_schedule_many = mocker.patch("gwf.plugins.run.Scheduler.schedule_many")

    args = ["-b", "testing", "run"]
    cli_runner.invoke(main, args)

    args, kwargs = mock_schedule_many.call_args
    assert len(args[0]) == 2
    assert {x.name for x in args[0]} == {"Target1", "Target2"}


def test_run_specified_target(cli_runner, mocker):
    mock_schedule_many = mocker.patch("gwf.plugins.run.Scheduler.schedule_many")

    args = ["-b", "testing", "run", "Target1"]
    cli_runner.invoke(main, args)

    args, kwargs = mock_schedule_many.call_args
    assert len(args[0]) == 1
    assert {x.name for x in args[0]} == {"Target1"}
