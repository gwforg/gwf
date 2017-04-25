from ..core import schedule_many
from ..cli import pass_backend, pass_graph
from ..exceptions import TargetDoesNotExistError

import click


@click.command()
@click.argument('targets', nargs=-1)
@pass_backend
@pass_graph
@click.pass_obj
def run(ctx, graph, backend, targets):
    """Run the specified workflow."""
    targets = []
    if not targets:
        targets = graph.endpoints()
    else:
        for name in targets:
            if name not in graph.targets:
                raise TargetDoesNotExistError(name)
            targets.append(graph.targets[name])
    schedule_many(graph, backend, targets)
