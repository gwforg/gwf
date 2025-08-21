from gwf.cli import main


def test_main_shows_usage_when_no_subcommand_is_given(cli_runner):
    result = cli_runner.invoke(main, [])
    assert result.output.startswith("Usage:")
