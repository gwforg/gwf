import json

from gwf.cli import main


def test_info_all_targets(cli_runner, simple_workflow):
    args = ["info"]
    result = cli_runner.invoke(main, args)
    doc = json.loads(result.output)

    assert "Target1" in doc
    assert "Target2" in doc
    assert "Target3" in doc


def test_info_single_target(cli_runner, simple_workflow):
    args = ["info", "Target1"]
    result = cli_runner.invoke(main, args)
    doc = json.loads(result.output)

    assert "Target1" in doc
    assert "Target2" not in doc
    assert "Target3" not in doc


def test_info_pretty_output(cli_runner, simple_workflow):
    args = ["info", "--format", "pretty"]
    result = cli_runner.invoke(main, args)

    assert "Target1" in result.output
    assert "Target2" in result.output
    assert "Target3" in result.output
    assert "Spec:" in result.output
