import logging
import os
from pkg_resources import iter_entry_points

import click
from click_plugins import with_plugins

from . import __version__
from .backends import Backend
from .conf import config
from .exceptions import ConfigurationError
from .utils import ColorFormatter, ensure_dir, get_latest_version

logger = logging.getLogger(__name__)


BASIC_FORMAT = "%(message)s"

ADVANCED_FORMAT = "%(levelname)s%(message)s"

LOGGING_FORMATS = {
    "warning": BASIC_FORMAT,
    "info": BASIC_FORMAT,
    "debug": ADVANCED_FORMAT,
    "error": BASIC_FORMAT,
}

VERBOSITY_LEVELS = ["warning", "debug", "info", "error"]


def get_level(level):
    return getattr(logging, level.upper())


def configure_logging(level_name):
    fmt = LOGGING_FORMATS[level_name]

    handler = logging.StreamHandler()
    handler.setFormatter(ColorFormatter(fmt=fmt))

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(get_level(level_name))


def _validate_choice(key, value, valid_values, human_values=None):
    if value not in valid_values:
        msg = 'Invalid value "{}" for key "{}", must be one of: {}.'
        values = valid_values if human_values is None else human_values
        raise ConfigurationError(msg.format(value, key, ", ".join(values)))


def _validate_bool(key, value):
    human_values = ("true", "yes", "false", "no")
    return _validate_choice(key, value, (True, False), human_values)


@config.validator("backend")
def validate_backend(value):
    return _validate_choice("backend", value, Backend.list())


@config.validator("verbose")
def validate_verbose(value):
    return _validate_choice("verbose", value, VERBOSITY_LEVELS)


@config.validator("no_color")
def validate_no_color(value):
    return _validate_bool("no_color", value)


@config.validator("check_updates")
def validate_check_updates(value):
    return _validate_bool("check_updates", value)


@with_plugins(iter_entry_points("gwf.plugins"))
@click.group(context_settings={"obj": {}})
@click.version_option(version=__version__)
@click.option("-f", "--file", default="workflow.py:gwf", help="Workflow/obj to load.")
@click.option(
    "-b",
    "--backend",
    type=click.Choice(Backend.list()),
    default=config["backend"],
    help="Backend used to run workflow.",
)
@click.option(
    "-v",
    "--verbose",
    type=click.Choice(VERBOSITY_LEVELS),
    default=config["verbose"],
    help="Verbosity level.",
)
@click.option(
    "--no-color/--use-color", default=None, help="Enable or disable output colors."
)
@click.pass_context
def main(ctx, file, backend, verbose, no_color):
    """A flexible, pragmatic workflow tool.

    See help for each command using the `--help` flag for that command:

        gwf status --help

    Shows help for the status command.
    """
    ensure_dir(".gwf")

    # If the --use-color/--no-color argument is not set, get a value from the
    # configuration file. If nothing has been configured, check if the NO_COLOR
    # environment variable has been set.
    if no_color is None:
        if config.get("no_color") is None:
            no_color = bool(os.getenv("NO_COLOR", False))
        else:
            no_color = config["no_color"]

    if no_color:
        # Hack for disabling all click colors. We basically lie to
        # click and pretend that the shell is never a TTY.
        click._compat.isatty = lambda s: False

    configure_logging(level_name=verbose)

    # Check for updates and show a warning if there's a new version.
    if config["check_updates"]:
        latest_version = get_latest_version()
        if latest_version is not None and latest_version != __version__:
            msg = (
                "You are currently running version {current_version}, but the "
                "latest version is {latest_version}. Consider updating. You "
                "can disable this check with `gwf config set check_updates no`."
            )
            logger.warning(
                msg.format(current_version=__version__, latest_version=latest_version)
            )

    ctx.obj = {"file": file, "backend": backend}
