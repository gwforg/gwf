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
    stdout_path, stderr_path = get_log_paths(ctx.working_dir, target)
    log_path = stderr_path if stderr else stdout_path
    try:
        # Targets may write output that is not valid UTF-8.
        log_file = log_path.open(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        stream = "standard error" if stderr else "standard output"
        raise click.ClickException(
            f"No {stream} log found for target {target}. "
            "Maybe the target has not run yet?"
        ) from None

    with log_file:
        if no_pager:
            for line in log_file:
                click.echo(line, nl=False)
        else:
            click.echo_via_pager(log_file)
