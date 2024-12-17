import logging
import os
import os.path
from contextlib import suppress

import click

from .. import Workflow
from ..backends import create_backend
from ..core import CachedFilesystem, Graph, get_spec_hashes, pass_context
from ..filtering import GroupFilter, NameFilter, filter_generic, filter_names
from ..scheduling import submit_workflow

logger = logging.getLogger(__name__)


def clean_logs(working_dir, graph):
    target_set = set(graph.targets.keys())
    log_files = set(
        os.path.splitext(os.path.basename(log_file))[0]
        for log_file in os.listdir(os.path.join(working_dir, ".gwf", "logs"))
    )

    for log_name in log_files.difference(target_set):
        logger.debug("Target %s does not exist, deleting", log_name)
        with suppress(OSError):
            os.remove(os.path.join(working_dir, ".gwf", "logs", log_name + ".stdout"))
            os.remove(os.path.join(working_dir, ".gwf", "logs", log_name + ".stderr"))


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-d", "--dry-run", is_flag=True, default=False)
@click.option("-f", "--force", is_flag=True, default=False)
@click.option("-n", "--no-deps", is_flag=True, default=False)
@click.option(
    "-g",
    "--group",
    multiple=True,
)
@pass_context
def run(ctx, targets, dry_run, force, no_deps, group):
    """Run the specified workflow."""
    workflow = Workflow.from_context(ctx)

    fs = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, fs)

    if ctx.config.get("clean_logs") and not dry_run:
        logger.debug("Cleaning unused log files...")
        clean_logs(ctx.working_dir, graph)

    with create_backend(
        ctx.backend, working_dir=ctx.working_dir, config=ctx.config
    ) as backend, get_spec_hashes(
        working_dir=ctx.working_dir, config=ctx.config
    ) as spec_hashes:
        filters = []
        if targets:
            filters.append(NameFilter(patterns=targets))
        if group:
            filters.append(GroupFilter(patterns=group))

        endpoints = set(filter_generic(targets=graph, filters=filters))

        # endpoints = filter_names(graph, targets) if targets else graph.endpoints()
        submit_workflow(
            endpoints,
            graph,
            fs,
            spec_hashes,
            backend,
            dry_run=dry_run,
            force=force,
            no_deps=no_deps,
        )
