import logging
import os

from ..core import graph_from_config
from ..filtering import NameFilter, EndpointFilter, filter_generic

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
@click.pass_obj
def clean(obj, targets, all):
    """Clean output files of targets.

    By default, only targets that are not endpoints will have their output files deleted. If you want to clean up output
    files from endpoints too, use the ``--all`` flag.
    """
    graph = graph_from_config(obj)

    filters = []
    if targets:
        filters.append(NameFilter(patterns=targets))
    if not all:
        filters.append(EndpointFilter(endpoints=graph.endpoints(), mode='exclude'))

    matches = filter_generic(targets=graph, filters=filters)
    for target in matches:
        logger.info('Deleting %s', target.name)
        for path in target.outputs:
            logging.info('Deleting output file "%s" from target "%s".', click.format_filename(path), target.name)
            _delete_file(path)
