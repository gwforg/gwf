from ..core import schedule_many
from ..cli import pass_backend, pass_graph
from ..exceptions import TargetDoesNotExistError

import click


@click.command()
@click.argument('targets', nargs=-1)
@click.option('-d', '--dry-run', is_flag=True, default=False)
@pass_backend
@pass_graph
def run(graph, backend, targets, dry_run):
    """Run the specified workflow."""
    targets = []
    if not targets:
        targets = graph.endpoints()
    else:
        for name in targets:
            if name not in graph.targets:
                raise TargetDoesNotExistError(name)
            targets.append(graph.targets[name])
    schedule_many(graph, backend, targets, dry_run=dry_run)
