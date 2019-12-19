import click
from ..core import Graph
from ..utils import touchfile


@click.command()
@click.pass_obj
def touch(obj):
    """Touch output files to update timestamps.

    Running this command touches all output files in the workflow such that
    their modification timestamp is updated. Touching is performed bottom-up
    such that, when done, all targets in the workflow will look completed.

    This is useful if one or more files were accidentially deleted, but you
    don't want to re-run the workflow to recreate them.
    """
    graph = Graph.from_config(obj)
    visited = set()
    for endpoint in graph.endpoints():
        for target in graph.dfs(endpoint):
            if target in visited:
                continue
            visited.add(target)
            for path in target.flattened_outputs():
                touchfile(path)
