import os.path

import pytest

from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt'])
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt'])
gwf.target('Target3', inputs=['a.txt'], outputs=['c.txt'])
"""


@pytest.fixture
def simple_workflow(tmpdir):
    workflow_file = tmpdir.join("workflow.py")
    workflow_file.write(SIMPLE_WORKFLOW)
    return tmpdir


@pytest.fixture(autouse=True)
def setup(simple_workflow):
    simple_workflow.join("a.txt").ensure()
    simple_workflow.join("b.txt").ensure()
    simple_workflow.join("c.txt").ensure()
    with simple_workflow.as_cwd():
        yield


def test_clean_output_from_non_endpoints(cli_runner):
    args = ["-b", "testing", "clean"]
    cli_runner.invoke(main, args, input="y\n")

    assert not os.path.exists("a.txt")
    assert os.path.exists("b.txt")
    assert os.path.exists("c.txt")


def test_clean_output_from_all_targets(cli_runner):
    args = ["-b", "testing", "clean", "--all"]
    cli_runner.invoke(main, args, input="y\n")

    assert not os.path.exists("a.txt")
    assert not os.path.exists("b.txt")
    assert not os.path.exists("c.txt")


def test_clean_output_from_single_endpoint_target(cli_runner):
    args = ["-b", "testing", "clean", "--all", "Target2"]
    cli_runner.invoke(main, args)

    assert os.path.exists("a.txt")
    assert not os.path.exists("b.txt")
    assert os.path.exists("c.txt")


def test_clean_output_from_two_targets(cli_runner):
    args = ["-b", "testing", "clean", "--all", "Target1", "Target2"]
    cli_runner.invoke(main, args)

    assert not os.path.exists("a.txt")
    assert not os.path.exists("b.txt")
