import click

from pathlib import Path

from ..backends import Backend
from ..conf import config_from_path


banner = r"""
  .-_'''-.   .--.      .--. ________
 '_( )_   \  |  |_     |  ||        |
|(_ o _)|  ' | _( )_   |  ||   .----'
. (_,_)/___| |(_ o _)  |  ||  _|____
|  |  .-----.| (_,_) \ |  ||_( )_   |
'  \  '-   .'|  |/    \|  |(_ o._)__|
 \  `-'`   | |  '  /\  `  ||(_,_)
  \        / |    /  \    ||   |
   `'-...-'  `---'    `---`'---'"""


workflow_template = '''from gwf import Workflow

from templates import *

gwf = Workflow()

gwf.target('ExampleTarget', inputs=[], outputs=[]) << """
echo hello world with gwf
"""
'''


@click.command()
def init():
    """Initialize a new gwf workflow.

    Running this command will create a new gwf workflow in an existing or new
    directory. This will generate a `workflow.py` file with an example
    workflow, as well as a `templates.py` file for your templates. It will also
    allow you to choose a backend.
    """
    width, height = click.get_terminal_size()
    max_width = max(len(line) for line in banner.splitlines())
    indent = int(width / 2 - max_width / 2)
    centered_banner = "\n".join((" " * indent) + line for line in banner.splitlines())

    click.echo()
    click.secho(centered_banner, fg="blue")
    click.echo()
    click.echo("Welcome!".center(width))
    click.echo()
    click.echo()

    project_dir = click.prompt(
        "Where do you want your project?",
        default=Path("."),
        value_proc=lambda x: Path(x),
    )
    backend = click.prompt(
        "Which backend do you want to use?",
        default="local",
        type=click.Choice(Backend.list()),
        show_choices=True,
    )

    click.echo("Generating project skeleton...")
    project_dir.mkdir(parents=True, exist_ok=True)
    project_dir.joinpath("workflow.py").write_text(workflow_template)
    project_dir.joinpath("templates.py").touch()

    click.echo("Setting project configuration...")

    config = config_from_path(project_dir.joinpath(".gwfconf.json"))
    config["backend"] = backend
    config.dump()

    click.echo()
    click.secho("Success!", fg="green")
    click.echo("Now go to {} and run 'gwf status' :-)".format(project_dir.resolve()))
