import click

from ..core import pass_context
from ..log_storage import get_log_paths


@click.command()
@click.argument("target")
@click.option("-e", "--stderr", is_flag=True)
@click.option("--no-pager", is_flag=True)
@pass_context
def logs(ctx, target, stderr, no_pager):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show
    standard error instead.
    """
    echo_func = click.echo if no_pager else click.echo_via_pager
    stdout_path, stderr_path = get_log_paths(ctx.working_dir, target)
    log_path = stderr_path if stderr else stdout_path
    echo_func(log_path.open(encoding="utf-8"))
