from functools import lru_cache
from pathlib import Path

import click

from .. import Workflow
from ..core import CachedFilesystem, Graph, get_spec_hashes, pass_context
from ..filtering import filter_names


def touch_workflow(endpoints, graph, spec_hashes):
    @lru_cache(maxsize=None)
    def _visit(target):
        for dep in graph.dependencies[target]:
            _visit(dep)

        spec_hashes.update(target)
        for path in target.flattened_outputs():
            Path(path).touch(exist_ok=True)

    for target in endpoints:
        _visit(target)


@click.command()
@click.argument("targets", nargs=-1)
@pass_context
def touch(ctx, targets):
    """Touch output files to update timestamps.

    Running this command touches all output files in the workflow such that
    their modification timestamp is updated. Touching is performed bottom-up
    such that, when done, all targets in the workflow will look completed. Spec
    hashes will also be "touched".

    This is useful if one or more files were accidentially deleted, but you
    don't want to re-run the workflow to recreate them.
    """
    workflow = Workflow.from_context(ctx)
    filesystem = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, filesystem)
    endpoints = filter_names(graph, targets) if targets else graph.endpoints()
    with get_spec_hashes(working_dir=ctx.working_dir, config=ctx.config) as spec_hashes:
        touch_workflow(endpoints, graph, spec_hashes)
