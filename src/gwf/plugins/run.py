import logging

import click

from .. import Workflow
from ..backends import create_backend
from ..log_storage import clean_logs
from ..core import CachedFilesystem, Graph, Status, get_spec_hashes, pass_context
from ..filtering import GroupFilter, NameFilter, StatusFilter, filter_generic
from ..scheduling import get_status_map, submit_workflow

logger = logging.getLogger(__name__)


_STATUS_NAMES = [s.name.lower() for s in Status]


def _status_names_to_enums(iterable):
    return [Status[s.upper()] for s in iterable]


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-d", "--dry-run", is_flag=True, default=False)
@click.option("-f", "--force", is_flag=True, default=False)
@click.option("-n", "--no-deps", is_flag=True, default=False)
@click.option("-g", "--group", multiple=True)
@click.option(
    "-s",
    "--status",
    type=click.Choice(_STATUS_NAMES),
    multiple=True,
)
@pass_context
def run(ctx, targets, dry_run, force, no_deps, group, status):
    """Run the specified workflow."""
    workflow = Workflow.from_context(ctx)

    fs = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets)

    if ctx.config.get("clean_logs") and not dry_run:
        logger.debug("Cleaning unused log files...")
        clean_logs(ctx.working_dir, graph)

    with create_backend(
        ctx.backend, working_dir=ctx.working_dir, config=ctx.config
    ) as backend:
        with get_spec_hashes(
            working_dir=ctx.working_dir, config=ctx.config
        ) as spec_hashes:
            target_states = get_status_map(
                graph,
                fs,
                spec_hashes,
                backend,
            )

            filters = []
            if targets:
                filters.append(NameFilter(patterns=targets))
            endpoints = set(filter_generic(targets=graph, filters=filters))

            filters = []
            if status:
                filters.append(
                    StatusFilter(
                        status_provider=target_states.get,
                        status=_status_names_to_enums(status),
                    )
                )
            if group:
                filters.append(GroupFilter(patterns=group))
            filtered_targets = set(filter_generic(targets=graph, filters=filters))
            filtered_graph = Graph.from_targets(filtered_targets)

            submit_workflow(
                endpoints,
                filtered_graph,
                fs,
                spec_hashes,
                backend,
                dry_run=dry_run,
                force=force,
                no_deps=no_deps,
            )
