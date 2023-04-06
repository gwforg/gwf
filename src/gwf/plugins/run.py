import logging
from contextlib import suppress

import click

from ..backends import Backend
from ..backends.exceptions import LogError
from ..conf import config
from ..core import CachedFilesystem, Graph, get_spec_hashes
from ..filtering import filter_names
from ..scheduling import submit_workflow
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
    fs = CachedFilesystem()
    workflow = Workflow.from_config(obj)
    graph = Graph.from_targets(workflow.targets, fs)

    backend_cls = Backend.from_config(obj)
    with backend_cls() as backend, get_spec_hashes(
        working_dir=workflow.working_dir, config=config
    ) as spec_hashes:
        if config.get("clean_logs") and not dry_run:
            logger.debug("Cleaning old log files...")
            clean_logs(graph, backend)

        endpoints = filter_names(graph, targets) if targets else graph.endpoints()
        submit_workflow(
            endpoints,
            graph,
            fs,
            spec_hashes,
            backend,
            dry_run=dry_run,
        )
