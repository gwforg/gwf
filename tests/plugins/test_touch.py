import os
import os.path

from gwf.cli import main


def test_touch_creates_files(cli_runner, local_backend, simple_workflow):
    cli_runner.invoke(main, ["touch"])

    assert os.path.exists("a.txt")
    assert os.path.exists("b.txt")
    assert os.path.exists("c.txt")
