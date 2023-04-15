import os.path

import pytest

from gwf.cli import main


@pytest.fixture(autouse=True)
def setup(simple_workflow):
    simple_workflow.join("a.txt").ensure()
    simple_workflow.join("b.txt").ensure()
    simple_workflow.join("c.txt").ensure()


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
