from functools import lru_cache

import click

from .. import Workflow
from ..core import CachedFilesystem, Graph, get_spec_hashes
from ..utils import touchfile


def touch_workflow(graph, spec_hashes):
    @lru_cache(maxsize=None)
    def _visit(target):
        for dep in graph.dependencies[target]:
            _visit(dep)

        spec_hashes.update(target)
        for path in target.flattened_outputs():
            touchfile(path)

    for target in graph.endpoints():
        _visit(target)


@click.command()
@click.pass_obj
def touch(obj):
    """Touch output files to update timestamps.

    Running this command touches all output files in the workflow such that
    their modification timestamp is updated. Touching is performed bottom-up
    such that, when done, all targets in the workflow will look completed. Spec
    hashes will also be "touched".

    This is useful if one or more files were accidentially deleted, but you
    don't want to re-run the workflow to recreate them.
    """
    workflow = Workflow.from_config(obj)
    filesystem = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, filesystem)
    with get_spec_hashes(
        working_dir=obj["working_dir"], config=obj["config"]
    ) as spec_hashes:
        touch_workflow(graph, spec_hashes)
