import logging
import os

from ..cli import pass_graph

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
@pass_graph
def clean(graph, targets, not_endpoints):
    """Clean output files of targets."""
    matched_targets = graph.find(targets) or graph.targets.values()

    if not_endpoints:
        matched_targets = [
            target
            for target in matched_targets
            if target not in graph.endpoints()
        ]

    for target in matched_targets:
        logger.info('Deleting %s', target.name)
        for path in target.outputs:
            logging.info('Deleting output file "%s" from target "%s".', click.format_filename(path), target.name)
            _delete_file(path)
