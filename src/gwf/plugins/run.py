import logging

import click

from ..backends import Backend, Status
from ..core import Graph, schedule
from ..filtering import filter_names
from ..workflow import Workflow

logger = logging.getLogger(__name__)


def submit(graph, scheduled, reasons, backend, dry_run):
    seen = set()
    for endpoint in graph.endpoints():
        for target in graph.dfs(endpoint):
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
                backend.submit_full(target, dependencies=scheduled[target])


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
        matched_targets = filter_names(graph, targets) if targets else graph.endpoints()
        subgraph = graph.subset(matched_targets)

        scheduled, reasons = schedule(matched_targets, subgraph)
        submit(subgraph, scheduled, reasons, backend, dry_run=dry_run)
