from gwf.cli import main


def test_set_get(cli_runner, simple_workflow):
    res = cli_runner.invoke(main, ["config", "set", "backend", "slurm"])
    assert res.exit_code == 0

    res = cli_runner.invoke(main, ["config", "get", "backend"])
    assert res.exit_code == 0
    assert res.output == "slurm\n"


def test_unset(cli_runner, simple_workflow):
    res = cli_runner.invoke(main, ["config", "unset", "backend"])
    assert res.exit_code == 0

    res = cli_runner.invoke(main, ["config", "get", "backend"])
    assert res.exit_code == 0
    assert res.output == "<not set>\n"
