import logging
import os
import os.path
import shutil
from pathlib import Path

import click
from click_plugins import with_plugins

from . import __version__
from .backends import guess_backend, list_backends
from .conf import VERBOSITY_LEVELS, FileConfig
from .exceptions import ConfigurationError
from .utils import ColorFormatter, entry_points, find_workflow

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


def _validate_choice(key, value, valid_values, human_values=None):
    if value not in valid_values:
        msg = 'Invalid value "{}" for key "{}", must be one of: {}.'
        values = valid_values if human_values is None else human_values
        raise ConfigurationError(msg.format(value, key, ", ".join(values)))


def _validate_bool(key, value):
    human_values = ("true", "yes", "false", "no")
    return _validate_choice(key, value, (True, False), human_values)


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
    type=click.Choice(VERBOSITY_LEVELS),
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
    configure_logging(level_name=verbose)

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
    working_dir.joinpath(".gwf").mkdir(exist_ok=True)
    working_dir.joinpath(".gwf", "logs").mkdir(exist_ok=True)

    config = FileConfig.load(working_dir.joinpath(".gwfconf.json"))

    @config.validator("backend")
    def validate_backend(value):
        return _validate_choice("backend", value, list_backends())

    @config.validator("verbose")
    def validate_verbose(value):
        return _validate_choice("verbose", value, VERBOSITY_LEVELS)

    @config.validator("no_color")
    def validate_no_color(value):
        return _validate_bool("no_color", value)

    @config.validator("use_spec_hashes")
    def validate_use_spec_hashes(value):
        return _validate_bool("use_spec_hashes", value)

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
        logger.debug("No backend was configured, guessing a backend instead")
        score, backend = guess_backend()
        logger.debug("Found backend '%s' with priority %s", backend, score)
    logger.debug("Using '%s' backend", backend)

    ctx.obj = {
        "config": config,
        "backend": backend,
        "working_dir": str(working_dir),
        "file": path,
        "obj_name": obj_name,
    }
