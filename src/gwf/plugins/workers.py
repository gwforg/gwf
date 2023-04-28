import multiprocessing
from threading import Thread

import click

from ..backends.local import Cluster


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
def workers(host, port, num_workers):
    """Start workers for the local backend."""
    cluster = Cluster(hostname=host, port=port, num_workers=num_workers)
    thread = Thread(target=cluster.start)
    try:
        thread.start()
        while thread.is_alive():
            thread.join(0.01)
    except KeyboardInterrupt:
        cluster.shutdown()
        thread.join()
