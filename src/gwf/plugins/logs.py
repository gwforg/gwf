import os.path

import click

from ..core import pass_context


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

    log_name = target + (".stderr" if stderr else ".stdout")
    log_path = os.path.join(ctx.logs_dir, log_name)
    echo_func(open(log_path).read())
