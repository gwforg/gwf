import os.path

import pytest

from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt']) << "echo hello world"
gwf.target('Target2', inputs=[], outputs=['b.txt']) << "echo world hello"
"""


@pytest.fixture
def simple_workflow(tmpdir):
    workflow_file = tmpdir.join('workflow.py')
    workflow_file.write(SIMPLE_WORKFLOW)
    return tmpdir


@pytest.fixture(autouse=True)
def setup(simple_workflow):
    simple_workflow.join('a.txt').ensure()
    simple_workflow.join('b.txt').ensure()
    with simple_workflow.as_cwd():
        yield


def test_clean_output_from_all_targets_by_default(cli_runner):
    args = ['-b', 'testing', 'clean']
    cli_runner.invoke(main, args)

    assert not os.path.exists('a.txt')
    assert not os.path.exists('b.txt')

def test_clean_output_from_single_target(cli_runner):
    args = ['-b', 'testing', 'clean', 'Target1']
    cli_runner.invoke(main, args)

    assert not os.path.exists('a.txt')
    assert os.path.exists('b.txt')

def test_clean_output_from_two_targets(cli_runner):
    args = ['-b', 'testing', 'clean', 'Target1', 'Target2']
    cli_runner.invoke(main, args)

    assert not os.path.exists('a.txt')
    assert not os.path.exists('b.txt')

def test_do_not_clean_outputs_from_endpoints(cli_runner):
    args = ['-b', 'testing', 'clean', '--not-endpoints', 'Target1', 'Target2']
    cli_runner.invoke(main, args)

    assert os.path.exists('a.txt')
    assert os.path.exists('b.txt')
