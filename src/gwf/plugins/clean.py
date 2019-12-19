import logging
import os
import os.path

from ..core import Graph
from ..filtering import NameFilter, EndpointFilter, filter_generic

import click

logger = logging.getLogger(__name__)


def _format_size(num, suffix="B"):
    # Implementation taken from:
    # https://stackoverflow.com/questions/1094841/reusable-library-to-get-human-readable-version-of-file-size
    for unit in ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi"]:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, "Yi", suffix)


def _delete_file(path):
    try:
        os.remove(path)
    except OSError:
        msg = "Error when attempting to delete {}, maybe the file did not exist"
        logger.debug(msg.format(path))


@click.command()
@click.argument("targets", nargs=-1)
@click.option("--all", is_flag=True, default=False)
@click.option(
    "-f", "--force", is_flag=True, default=False, help="Do not ask for confirmation."
)
@click.pass_obj
def clean(obj, targets, all, force):
    """Clean output files of targets.

    By default, only targets that are not endpoints will have their output files
    deleted. If you want to clean up output files from endpoints too, use the
    ``--all`` flag.
    """
    graph = Graph.from_config(obj)

    filters = []
    if targets:
        filters.append(NameFilter(patterns=targets))
    if not all:
        filters.append(EndpointFilter(endpoints=graph.endpoints(), mode="exclude"))

    matches = list(filter_generic(targets=graph, filters=filters))

    total_size = sum(
        os.path.getsize(path)
        if os.path.exists(path) and path not in target.protected
        else 0
        for target in matches
        for path in target.flattened_outputs()
    )

    logger.info("Will delete %s of files!", _format_size(total_size))

    if not targets and not force:
        click.confirm(
            "This will delete all unprotected output files from non-endpoint targets! Do you want to continue?",
            abort=True,
        )

    for target in matches:
        logger.info("Deleting output files of %s", target.name)
        for path in target.flattened_outputs():
            if path in target.protected:
                logging.debug(
                    'Skpping deleting file "%s" from target "%s" because it is protected',
                    click.format_filename(path),
                    target.name,
                )
                continue

            logging.info(
                'Deleting file "%s" from target "%s"',
                click.format_filename(path),
                target.name,
            )
            _delete_file(path)
