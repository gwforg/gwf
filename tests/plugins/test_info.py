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
