import multiprocessing

import click

from ..config import conf
from ..backends.local import Server


@click.command()
@click.option(
    '-n',
    '--num-workers',
    type=int,
    default=conf.get('local.num_workers', multiprocessing.cpu_count()),
    help='Number of workers to spawn.'
)
@click.option(
    '-p',
    '--port',
    type=int,
    default=conf.get('local.port', 12345),
    help='Port that workers will listen on.'
)
@click.option(
    '-h',
    '--host',
    default=conf.get('local.host', 'localhost'),
    help='Host that workers will bind to.'
)
def workers(host, port, num_workers):
    """Start workers for the local backend."""
    server = Server(hostname=host, port=port, num_workers=num_workers)
    server.start()
