import logging
from contextlib import suppress

import click

from ..backends import Backend, Status
from ..backends.exceptions import LogError
from ..conf import config
from ..core import Graph, dump_spec_hashes, load_spec_hashes
from ..filtering import filter_names
from ..scheduling import linearize_plan, schedule_workflow
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


def submit_workflow(matched_targets, reasons, backend, dry_run):
    seen = set()
    for endpoint in matched_targets:
        reason = reasons[endpoint]
        plan = linearize_plan(reason)
        for reason, deps in plan:
            target = reason.target
            if target in seen:
                continue

            if backend.status(target) != Status.UNKNOWN:
                logger.debug("Target %s already submitted", target.name)
                continue

            if dry_run:
                logger.info(
                    "Would submit target %s (%s)",
                    target,
                    reason,
                )
            else:
                logger.info("Submitting target %s (%s)", target.name, reason.reason())
                backend.submit_full(target, dependencies=deps)
            seen.add(target)


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-d", "--dry-run", is_flag=True, default=False)
@click.pass_obj
def run(obj, targets, dry_run):
    """Run the specified workflow."""
    workflow = Workflow.from_config(obj)
    graph = Graph.from_targets(workflow.targets)

    matched_targets = filter_names(graph, targets) if targets else graph.endpoints()

    spec_hashes = None
    if config.get("use_spec_hashes"):
        spec_hashes = load_spec_hashes(
            workflow.working_dir,
            graph,
        )

    reasons = schedule_workflow(
        graph,
        spec_hashes=spec_hashes,
        endpoints=matched_targets,
    )

    backend_cls = Backend.from_config(obj)
    with backend_cls() as backend:
        if config.get("clean_logs") and not dry_run:
            logger.debug("Cleaning old log files...")
            clean_logs(graph, backend)

        submit_workflow(matched_targets, reasons, backend, dry_run)

    if spec_hashes is not None and not dry_run:
        dump_spec_hashes(workflow.working_dir, graph)
