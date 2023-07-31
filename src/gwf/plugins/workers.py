import multiprocessing

import click

from ..backends.local import start_cluster
from ..core import pass_context


@click.command()
@click.option(
    "-n",
    "--num-workers",
    type=int,
    default=multiprocessing.cpu_count(),
    help="Number of workers to spawn.",
)
@click.option(
    "-p",
    "--port",
    type=int,
    default=12345,
    help="Port that workers will listen on.",
)
@click.option(
    "-h",
    "--host",
    default="localhost",
    help="Host that workers will bind to.",
)
@pass_context
def workers(ctx, host, port, num_workers):
    """Start workers for the local backend."""
    start_cluster(ctx.working_dir, num_workers, host, port)
