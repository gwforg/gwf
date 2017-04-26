import multiprocessing

import click

from ..backends.local import Server


@click.command()
@click.option('-n', '--num-workers', type=int, default=multiprocessing.cpu_count(), help='Number of workers to spawn.')
@click.option('-p', '--port', type=int, default=12345, help='Port that workers will listen on.')
@click.pass_context
def workers(ctx, port, num_workers):
    """Start workers for the local backend."""
    server = Server(port=port, num_workers=num_workers)
    server.start()
