import io
import logging

import click
import pytest

from gwf.cli import ADVANCED_FORMAT, main, set_colors
from gwf.utils import ColorFormatter


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


class _FakeTerminal(io.StringIO):
    def isatty(self):
        return True


@pytest.mark.parametrize(
    ("stream", "enabled", "ctx_color", "log_color"),
    [
        (_FakeTerminal(), True, None, True),
        (_FakeTerminal(), False, False, False),
        (io.StringIO(), True, None, False),
    ],
)
def test_set_colors(stream, enabled, ctx_color, log_color):
    ctx = click.Context(main)
    handler = logging.StreamHandler(stream)
    handler.setFormatter(ColorFormatter(fmt=ADVANCED_FORMAT))

    set_colors(ctx, handler, enabled)

    assert ctx.color is ctx_color
    assert handler.formatter.color is log_color


@pytest.mark.parametrize("color", [True, False])
def test_color_formatter(color):
    formatter = ColorFormatter(fmt=ADVANCED_FORMAT, color=color)
    record = logging.LogRecord("gwf", logging.DEBUG, __file__, 1, "hello", (), None)

    output = formatter.format(record)

    assert ("\x1b[" in output) is color
    assert click.unstyle(output) == "debug     hello"


@pytest.mark.parametrize(("args", "colored"), [([], True), (["--no-color"], False)])
def test_no_color_disables_colors_in_command_output(
    cli_runner, simple_workflow, args, colored
):
    # color=True makes the test runner behave like a terminal.
    result = cli_runner.invoke(main, [*args, "info", "-f", "pretty"], color=True)

    assert result.exit_code == 0, result.output
    assert ("\x1b[" in result.output) is colored
