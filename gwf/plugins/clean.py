import logging
import os

from ..cli import pass_graph
from ..filtering import Criteria, NameFilter, EndpointFilter, filter

import click

logger = logging.getLogger(__name__)


def _delete_file(path):
    try:
        os.remove(path)
    except OSError:
        pass


@click.command()
@click.argument('targets', nargs=-1)
@click.option('--all', is_flag=True, default=False)
@pass_graph
def clean(graph, targets, all):
    """Clean output files of targets.

    By default, only targets that are not endpoints will have their output files deleted. If you want to clean up output
    files from endpoints too, use the ``--all`` flag.
    """
    matches = filter(
        targets=graph.targets.values(),
        criteria=Criteria(targets=targets, all=all),
        filters=[
            NameFilter(),
            EndpointFilter(endpoints=graph.endpoints(), negate=True),
        ]
    )

    for target in matches:
        logger.info('Deleting %s', target.name)
        for path in target.outputs:
            logging.info('Deleting output file "%s" from target "%s".', click.format_filename(path), target.name)
            _delete_file(path)
