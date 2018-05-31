import pytest

from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow
gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << "echo hello world"
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt']) << "echo world hello"
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


def test_status_shows_all_targets(cli_runner):
    result = cli_runner.invoke(main, ["-b", "testing", "status"])
    assert "Target2" in result.output
    assert "Target1" in result.output


def test_status_shows_one_named_target(cli_runner):
    result = cli_runner.invoke(main, ["-b", "testing", "status", "Target1"])
    assert "Target2" not in result.output
    assert "Target1" in result.output


def test_status_shows_two_named_targets(cli_runner):
    result = cli_runner.invoke(main, ["-b", "testing", "status", "Target1", "Target2"])
    assert "Target2" in result.output
    assert "Target1" in result.output


def test_status_shows_only_endpoint_targets(cli_runner):
    result = cli_runner.invoke(main, ["-b", "testing", "status", "--endpoints"])
    assert "Target2" in result.output
    assert "Target1" not in result.output
