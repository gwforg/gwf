import logging
import os
import shutil
from importlib.metadata import entry_points
from pathlib import Path

import click
from click_plugins import with_plugins


from . import __version__
from .backends import guess_backend, list_backends
from .conf import FileConfig
from .core import Context
from .log_storage import init_log_storage
from .utils import ColorFormatter, find_workflow

logger = logging.getLogger(__name__)


BASIC_FORMAT = "%(message)s"

ADVANCED_FORMAT = "%(levelname)s%(message)s"

LOGGING_FORMATS = {
    "warning": BASIC_FORMAT,
    "info": BASIC_FORMAT,
    "debug": ADVANCED_FORMAT,
    "error": BASIC_FORMAT,
}

BANNER = r"""
  .-_'''-.   .--.      .--. ________
 '_( )_   \  |  |_     |  ||        |
|(_ o _)|  ' | _( )_   |  ||   .----'
. (_,_)/___| |(_ o _)  |  ||  _|____
|  |  .-----.| (_,_) \ |  ||_( )_   |
'  \  '-   .'|  |/    \|  |(_ o._)__|
 \  `-'`   | |  '  /\  `  ||(_,_)
  \        / |    /  \    ||   |
   `'-...-'  `---'    `---`'---'"""


WORKFLOW_TEMPLATE = '''from gwf import Workflow

from templates import *

gwf = Workflow()

gwf.target('ExampleTarget', inputs=[], outputs=[]) << """
echo hello world with gwf
"""
'''


def get_level(level):
    return getattr(logging, level.upper())


def configure_logging(level_name):
    fmt = LOGGING_FORMATS[level_name]

    handler = logging.StreamHandler()
    handler.setFormatter(ColorFormatter(fmt=fmt))

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(get_level(level_name))
    return handler


def set_colors(ctx, log_handler, enabled):
    """Enable or disable colors in command output and log messages."""
    # Subcommand contexts inherit the color setting, so this applies to all
    # output from click.echo() and friends. None lets click detect a terminal.
    ctx.color = None if enabled else False
    # Log messages are not written through click, so detect a terminal here.
    log_handler.formatter.color = enabled and log_handler.stream.isatty()


def init(project_dir):
    """Initialize a new gwf workflow.

    Running this command will create a new gwf workflow in an existing or new
    directory. This will generate a `workflow.py` file with an example
    workflow, as well as a `templates.py` file for your templates. It will also
    allow you to choose a backend.
    """
    width, height = shutil.get_terminal_size()
    max_width = max(len(line) for line in BANNER.splitlines())
    indent = int(width / 2 - max_width / 2)
    centered_banner = "\n".join((" " * indent) + line for line in BANNER.splitlines())

    click.echo()
    click.secho(centered_banner, fg="blue")
    click.echo()
    click.echo("Welcome!".center(width))
    click.echo()
    click.echo()

    _, backend_guess = guess_backend()
    backend = click.prompt(
        "Which backend do you want to use?",
        default=backend_guess,
        type=click.Choice(list_backends()),
        show_choices=True,
    )

    click.echo("Generating project skeleton...")
    project_dir.mkdir(parents=True, exist_ok=True)
    project_dir.joinpath("workflow.py").write_text(WORKFLOW_TEMPLATE)
    project_dir.joinpath("templates.py").touch()

    click.echo("Setting project configuration...")

    config = FileConfig.load(project_dir.joinpath(".gwfconf.json"))
    config["backend"] = backend
    config.dump()

    click.echo()
    click.secho("Success!", fg="green")
    click.echo("Now go to {} and run 'gwf status' :-)".format(project_dir.resolve()))


@with_plugins(entry_points(group="gwf.plugins"))
@click.group(context_settings={"obj": {}})
@click.version_option(version=__version__)
@click.option("-f", "--file", default="workflow.py:gwf", help="Workflow/obj to load.")
@click.option(
    "-b",
    "--backend",
    type=click.Choice(list_backends()),
    help="Backend used to run workflow.",
)
@click.option(
    "-v",
    "--verbose",
    type=click.Choice(["warning", "debug", "info", "error"]),
    default="info",
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
    log_handler = configure_logging(level_name=verbose)

    # Colors can be disabled with --no-color, in the configuration file, or
    # with the NO_COLOR environment variable, in that order of precedence.
    # Decide now so that early log messages are covered, and apply the
    # configuration file once it has been loaded.
    if no_color is None:
        set_colors(ctx, log_handler, enabled=not os.getenv("NO_COLOR"))
    else:
        set_colors(ctx, log_handler, enabled=not no_color)

    try:
        path, obj_name = find_workflow(file)
        working_dir = path.parent
    except FileNotFoundError:
        click.confirm(
            "Could not find a workflow file! Do you want to create one in this directory?",
            abort=True,
        )
        working_dir = Path.cwd()
        path, obj_name = working_dir.joinpath("workflow.py"), "gwf"
        init(working_dir)

    # Instantiate workflow config directory.
    config_dir = working_dir.joinpath(".gwf")
    config_dir.mkdir(exist_ok=True)

    # The configuration lives next to the workflow file, where `gwf init` and
    # all released versions of gwf put it.
    config = FileConfig.load(working_dir.joinpath(".gwfconf.json"))

    init_log_storage(working_dir)

    if no_color is None and config.get("no_color") is not None:
        set_colors(ctx, log_handler, enabled=not config["no_color"])

    backend = backend or config.get("backend")
    if backend is None:
        logger.debug("No backend was configured, guessing a backend instead")
        score, backend = guess_backend()
        logger.debug("Found backend '%s' with priority %s", backend, score)
    logger.debug("Using '%s' backend", backend)

    ctx.obj = Context(
        working_dir=str(working_dir),
        config=config,
        backend=backend,
        workflow_file=path,
        workflow_obj=obj_name,
    )
