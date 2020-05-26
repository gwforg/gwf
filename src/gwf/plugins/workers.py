import multiprocessing

import click

from ..conf import config
from ..backends.logmanager import FileLogManager
from ..backends.local import start_server


@click.command()
@click.option(
    "-n",
    "--max-cores",
    type=int,
    default=config.get("local.max_cores", multiprocessing.cpu_count()),
    help="Maximum number of cores to allocate.",
)
@click.option(
    "-p",
    "--port",
    type=int,
    default=config.get("local.port", 12345),
    help="Port that workers will listen on.",
)
@click.option(
    "-h",
    "--host",
    default=config.get("local.host", "localhost"),
    help="Host that workers will bind to.",
)
def workers(host, port, max_cores):
    """Start workers for the local backend."""
    log_manager = FileLogManager()
    start_server(log_manager, host=host, port=port, max_cores=max_cores)
