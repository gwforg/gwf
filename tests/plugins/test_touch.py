import os
import os.path

import pytest

from gwf.cli import main


SIMPLE_WORKFLOW = """from gwf import Workflow

gwf = Workflow()
gwf.target('Target1', inputs=[], outputs=['a.txt'])
gwf.target('Target2', inputs=['a.txt'], outputs=['b.txt'])
gwf.target('Target3', inputs=['b.txt'], outputs=['c.txt'])
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


def test_touch_creates_files(cli_runner):
    cli_runner.invoke(main, ["touch"])

    assert os.path.exists("a.txt")
    assert os.path.exists("b.txt")
    assert os.path.exists("c.txt")

    stat_a = os.stat("a.txt").st_ctime
    stat_b = os.stat("b.txt").st_ctime
    stat_c = os.stat("c.txt").st_ctime

    assert stat_a <= stat_b
    assert stat_b <= stat_c
