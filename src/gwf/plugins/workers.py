import multiprocessing

import click

from ..backends.local import start_cluster
from ..core import pass_context


@click.command()
@click.option(
    "-n",
    "--max-cores",
    "max_cores",
    type=click.IntRange(min=1),
    default=multiprocessing.cpu_count(),
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
