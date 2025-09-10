import os
import os.path

from gwf.cli import main


def test_touch_creates_files_with_no_create_missing(
    cli_runner, local_backend, simple_workflow
):
    cli_runner.invoke(main, ["touch"])

    assert not os.path.exists("a.txt")
    assert not os.path.exists("b.txt")
    assert not os.path.exists("c.txt")


def test_touch_creates_files_with_create_missing(
    cli_runner, local_backend, simple_workflow
):
    cli_runner.invoke(main, ["touch", "--create-missing"])

    assert os.path.exists("a.txt")
    assert os.path.exists("b.txt")
    assert os.path.exists("c.txt")
