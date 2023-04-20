import logging
import os
import os.path

import click
from click_plugins import with_plugins

from . import Workflow, __version__
from .backends import Backend
from .conf import VERBOSITY_LEVELS, config
from .utils import ColorFormatter, ensure_dir, entry_points, find_workflow

logger = logging.getLogger(__name__)


BASIC_FORMAT = "%(message)s"

ADVANCED_FORMAT = "%(levelname)s%(message)s"

LOGGING_FORMATS = {
    "warning": BASIC_FORMAT,
    "info": BASIC_FORMAT,
    "debug": ADVANCED_FORMAT,
    "error": BASIC_FORMAT,
}


def get_level(level):
    return getattr(logging, level.upper())


def configure_logging(level_name):
    fmt = LOGGING_FORMATS[level_name]

    handler = logging.StreamHandler()
    handler.setFormatter(ColorFormatter(fmt=fmt))

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(get_level(level_name))
    return root


@with_plugins(entry_points(group="gwf.plugins"))
@click.group(context_settings={"obj": {}})
@click.version_option(version=__version__)
@click.option("-f", "--file", default="workflow.py:gwf", help="Workflow/obj to load.")
@click.option(
    "-b",
    "--backend",
    type=click.Choice(Backend.list()),
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
    configure_logging(level_name=verbose)

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

    if backend is None:
        backend = config.get("backend")

    if backend is None:
        logging.debug("No backend was configured, guessing a backend instead")
        score, backend = Backend.guess()
        logger.debug("Found backend '%s' with priority %s", backend, score)
    logger.debug("Using '%s' backend", backend)

    path, obj_name = find_workflow(file)

    # Instantiate workflow config directory.
    path.parent.joinpath(".gwf").mkdir(exist_ok=True)
    path.parent.joinpath(".gwf", "logs").mkdir(exist_ok=True)

    ctx.obj = {
        "backend": backend,
        "working_dir": path.parent,
        "file": path,
        "obj_name": obj_name,
    }
