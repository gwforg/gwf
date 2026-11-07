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


def test_unset_key_with_only_a_default(cli_runner, simple_workflow):
    res = cli_runner.invoke(main, ["config", "unset", "clean_logs"])
    assert res.exit_code == 0, res.output

    res = cli_runner.invoke(main, ["config", "get", "clean_logs"])
    assert res.exit_code == 0
    assert res.output == "True\n"
