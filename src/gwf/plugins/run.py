import logging
from contextlib import suppress

import click

from ..backends import Backend, Status
from ..backends.exceptions import LogError
from ..core import Graph, schedule
from ..filtering import filter_names
from ..workflow import Workflow

logger = logging.getLogger(__name__)


def clean_logs(graph, backend):
    target_set = set(graph.targets.keys())
    log_files = set(backend.log_manager.list())

    for target_name in log_files.difference(target_set):
        logger.debug("Target %s does not exist, deleting log files...", target_name)

        with suppress(LogError):
            backend.log_manager.remove_stdout(target_name)
            backend.log_manager.remove_stderr(target_name)


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-d", "--dry-run", is_flag=True, default=False)
@click.pass_obj
def run(obj, targets, dry_run):
    """Run the specified workflow."""
    workflow = Workflow.from_config(obj)
    graph = Graph.from_targets(workflow.targets)
    backend_cls = Backend.from_config(obj)

    with backend_cls() as backend:
        if not dry_run:
            logger.debug("Cleaning old log files...")
            clean_logs(graph, backend)

        matched_targets = filter_names(graph, targets) if targets else graph.endpoints()
        subgraph = graph.subset(matched_targets)

        scheduled, reasons = schedule(matched_targets, subgraph)

        seen = set()
        for endpoint in subgraph.endpoints():
            for target in subgraph.dfs(endpoint):
                if target in seen:
                    continue
                seen.add(target)

                if target not in scheduled:
                    logger.debug(reasons[target])
                    continue

                if backend.status(target) != Status.UNKNOWN:
                    logger.debug("Target %s already submitted", target.name)
                    continue

                logger.debug("Target %s", reasons[target])
                if dry_run:
                    logger.info("Would submit target %s", target.name)
                else:
                    logger.info("Submitting target %s", target.name)
                    backend.submit(target, dependencies=scheduled[target])
