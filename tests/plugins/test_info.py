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


STRUCTURED_WORKFLOW = """from pathlib import Path
from gwf import Workflow

gwf = Workflow()
gwf.target(
    'Foo',
    inputs={'reads': [Path('r1.fq'), Path('r2.fq')], 'ref': 'ref.fa'},
    outputs='out.bam',
) << 'touch out.bam'
"""


def test_info_json_with_path_objects_and_named_inputs(cli_runner, tmpdir):
    tmpdir.join("workflow.py").write(STRUCTURED_WORKFLOW)

    with tmpdir.as_cwd():
        result = cli_runner.invoke(main, ["info"])

    assert result.exit_code == 0, result.output
    doc = json.loads(result.output)
    assert doc["Foo"]["inputs"] == {"reads": ["r1.fq", "r2.fq"], "ref": "ref.fa"}
    assert doc["Foo"]["outputs"] == "out.bam"


def test_info_pretty_with_path_objects_and_named_inputs(cli_runner, tmpdir):
    tmpdir.join("workflow.py").write(STRUCTURED_WORKFLOW)

    with tmpdir.as_cwd():
        result = cli_runner.invoke(main, ["info", "--format", "pretty"])

    assert result.exit_code == 0, result.output
    lines = [line.strip() for line in result.output.splitlines()]
    inputs = lines[lines.index("Inputs:") + 1 : lines.index("Outputs:")]
    outputs = lines[lines.index("Outputs:") + 1 : lines.index("Dependents:")]
    assert sorted(inputs) == ["r1.fq", "r2.fq", "ref.fa"]
    assert outputs == ["out.bam"]
