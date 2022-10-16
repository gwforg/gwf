import logging
import os.path
from contextlib import suppress

import click

from ..backends import Backend, Status
from ..backends.exceptions import LogError
from ..conf import config
from ..core import Graph, Reason, hash_spec, schedule
from ..filtering import filter_names
from ..utils import PersistableDict
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
        if config.get("clean_logs") and not dry_run:
            logger.debug("Cleaning old log files...")
            clean_logs(graph, backend)

        matched_targets = filter_names(graph, targets) if targets else graph.endpoints()
        subgraph = graph.subset(matched_targets)

        spec_hashes = None
        if config.get("use_spec_hashing"):
            spec_hashes = PersistableDict(os.path.join(".gwf", "spec-hashes.json"))

        plans = schedule(matched_targets, spec_hashes=spec_hashes, graph=subgraph)

        seen = set()

        def submit(reason):
            if reason.target not in seen and reason.scheduled:
                dependencies = set()
                if isinstance(reason, Reason.DependencyScheduled):
                    for dep in reason.dependencies:
                        submit(dep)
                        dependencies.add(dep.target)

                target = reason.target
                if backend.status(target) != Status.UNKNOWN:
                    logger.debug("Target %s already submitted", target.name)
                    return

                if dry_run:
                    logger.info(
                        "Would submit target %s (%s)",
                        target,
                        reason,
                    )
                else:
                    logger.info(
                        "Submitting target %s (%s)", target.name, reason.reason()
                    )
                    backend.submit_full(target, dependencies=dependencies)

        for _, plan in plans.items():
            submit(plan)

        if spec_hashes is not None and not dry_run:
            for target in graph.targets.values():
                spec_hashes[target.name] = hash_spec(target.spec)
            spec_hashes.persist()
