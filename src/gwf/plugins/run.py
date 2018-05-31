import logging

from ..backends import backend_from_config
from ..core import Scheduler, graph_from_config
from ..filtering import filter_names

import click


logger = logging.getLogger(__name__)


def clean_logs(graph, backend):
    target_set = set(graph.targets.keys())
    log_files = set(backend.log_manager.list())

    for target_name in log_files:
        if target_name not in target_set:
            logger.debug("Target %s does not exist, deleting log files...", target_name)

            backend.log_manager.remove_stdout(target_name)
            backend.log_manager.remove_stderr(target_name)


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-d", "--dry-run", is_flag=True, default=False)
@click.pass_obj
def run(obj, targets, dry_run):
    """Run the specified workflow."""
    graph = graph_from_config(obj)

    backend_cls = backend_from_config(obj)
    with backend_cls() as backend:
        if not dry_run:
            logger.debug("Cleaning old log files...")
            clean_logs(graph, backend)

        if targets:
            matched_targets = filter_names(graph, targets)
        else:
            matched_targets = list(graph)

        scheduler = Scheduler(graph=graph, backend=backend, dry_run=dry_run)
        scheduler.schedule_many(matched_targets)
