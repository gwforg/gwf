from gwf.cli import main


def test_main_shows_usage_when_no_subcommand_is_given(cli_runner):
    result = cli_runner.invoke(main, [])
    assert result.output.startswith("Usage:")


def test_init_configures_backend_used_by_later_commands(cli_runner, tmpdir):
    with tmpdir.as_cwd():
        # Create a project via the init prompt, choosing a backend that would
        # not be guessed here, then read the setting in the same invocation.
        result = cli_runner.invoke(
            main, ["config", "get", "backend"], input="y\nslurm\n"
        )

    assert result.exit_code == 0, result.output
    assert result.output.splitlines()[-1] == "slurm"


def test_config_next_to_workflow_file_is_used(cli_runner, simple_workflow):
    simple_workflow.join(".gwfconf.json").write('{"backend": "slurm"}\n')

    result = cli_runner.invoke(main, ["config", "get", "backend"])

    assert result.exit_code == 0, result.output
    assert result.output == "slurm\n"
