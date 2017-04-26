import logging
import os

from ..cli import pass_backend, pass_graph
from ..exceptions import TargetDoesNotExistError

import click

logger = logging.getLogger(__name__)


def _delete_file(path):
    try:
        os.remove(path)
    except OSError:
        pass


@click.command()
@click.argument('targets', nargs=-1)
@click.option('--not-endpoints', is_flag=True, default=False)
@pass_backend
@pass_graph
@click.pass_context
def clean(ctx, graph, backend, targets, not_endpoints):
    """Clean output files of targets."""
    targets = list(targets)
    if not targets:
        targets.extend(graph.targets.values())
    else:
        for name in targets:
            if name not in graph.targets:
                raise TargetDoesNotExistError(name)
            targets.append(graph.targets[name])

    if not_endpoints:
        for endpoint in graph.endpoints():
            targets.remove(endpoint)

    for target in targets:
        logger.info('Deleting %s', target.name)
        for path in target.outputs:
            logging.info('Deleting output file "%s" from target "%s".', click.format_filename(path), target.name)
            _delete_file(path)

