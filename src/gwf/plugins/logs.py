import click

from .. import Workflow
from ..backends import create_backend
from ..exceptions import GWFError


@click.command()
@click.argument("target")
@click.option("-e", "--stderr", is_flag=True)
@click.option("--no-pager", is_flag=True)
@click.pass_obj
def logs(obj, target, stderr, no_pager):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show
    standard error instead.
    """
    workflow = Workflow.from_config(obj)

    try:
        target = workflow.targets[target]
    except KeyError as exc:
        raise GWFError(f"Target {target} not found in the workflow.") from exc

    with create_backend(
        obj["backend"], working_dir=obj["working_dir"], config=obj["config"]
    ) as backend:
        log_file = backend.logs(target, stderr=stderr)
        log_contents = log_file.read()
        log_file.close()

    echo_func = click.echo if no_pager else click.echo_via_pager
    echo_func(log_contents)
