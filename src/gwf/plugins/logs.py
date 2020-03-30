import click

from ..backends import Backend
from ..workflow import Workflow
from ..exceptions import WorkflowError


@click.command()
@click.argument("target")
@click.option("-e", "--stderr", is_flag=True)
@click.option("--no-pager", is_flag=True)
@click.pass_obj
def logs(obj, target, stderr, no_pager):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show standard error instead.
    """
    workflow = Workflow.from_config(obj)
    backend_cls = Backend.from_config(obj)

    if target not in workflow.targets:
        raise WorkflowError('Target "{}" is not found in the workflow.'.format(target))

    log_file = backend_cls.logs(workflow.targets[target], stderr=stderr)
    log_contents = log_file.read()
    log_file.close()

    echo_func = click.echo if no_pager else click.echo_via_pager
    echo_func(log_contents)
