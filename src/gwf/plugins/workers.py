import click

from ..backends.local import get_total_cores, start_cluster
from ..core import pass_context


@click.command()
@click.option(
    "-n",
    "--max-cores",
    "max_cores",
    type=click.IntRange(min=1),
    # Evaluated when the command runs, not when gwf is imported.
    default=get_total_cores,
    help="Maximum number of cores to use.",
)
@click.option(
    "-m",
    "--max-memory",
    "max_memory",
    help="Maximum memory to use, for example 16g. Defaults to system memory.",
)
@pass_context
def workers(ctx, max_cores, max_memory):
    """Start the local scheduler."""
    start_cluster(ctx.working_dir, max_cores, max_memory=max_memory)
