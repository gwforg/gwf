import pytest

from gwf.cli import main
from gwf.backends.local import Client, LocalBackend, LocalStatus

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


def test_run_all_targets(cli_runner, local_backend):
    result = cli_runner.invoke(main, ["-v", "debug", "run"])
    assert result.exit_code == 0

    with LocalBackend() as backend:
        target1_id = backend.get_task_id("Target1")
        target2_id = backend.get_task_id("Target2")

    client = Client()
    client.connect()
    assert client.status(target1_id) == LocalStatus.SUBMITTED
    assert client.status(target2_id) == LocalStatus.SUBMITTED
    client.close()


def test_run_specified_target(cli_runner, local_backend):
    result = cli_runner.invoke(main, ["run", "Target1"])
    assert result.exit_code == 0

    with LocalBackend() as backend:
        task_id = backend.get_task_id("Target1")

    client = Client()
    client.connect()
    assert client.status(task_id) == LocalStatus.SUBMITTED
    client.close()
