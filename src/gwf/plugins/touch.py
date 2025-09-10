from functools import lru_cache
from pathlib import Path

import click

from .. import Workflow
from ..core import CachedFilesystem, Graph, get_spec_hashes, pass_context
from ..filtering import filter_names


def touch_workflow(endpoints, graph, spec_hashes, create_missing):
    @lru_cache(maxsize=None)
    def _visit(target):
        for dep in graph.dependencies[target]:
            _visit(dep)

        spec_hashes.update(target)
        for path in target.flattened_outputs():
            path = Path(path)
            if path.exists() or create_missing:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.touch(exist_ok=True)

    for target in endpoints:
        _visit(target)


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-c", "--create-missing", is_flag=True, default=False)
@pass_context
def touch(ctx, targets, create_missing):
    """Touch output files to update timestamps.

    Running this command touches all, existing output files in the workflow such
    that their modification timestamp is updated. Touching is performed
    bottom-up such that, when done, all targets in the workflow will look
    completed. Spec hashes will also be "touched".

    This is useful if one or more files were accidentially deleted, but you
    don't want to re-run the workflow to re-create them.

    By default, only files that already exist will be touched/updated. If `-c`
    is given, missing files will also be created.
    """
    workflow = Workflow.from_context(ctx)
    filesystem = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, filesystem)
    endpoints = filter_names(graph, targets) if targets else graph.endpoints()
    with get_spec_hashes(working_dir=ctx.working_dir, config=ctx.config) as spec_hashes:
        touch_workflow(endpoints, graph, spec_hashes, create_missing)
