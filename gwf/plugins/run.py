from ..backends import backend_from_config
from ..core import Scheduler, graph_from_config
from ..filtering import filter_names

import click


@click.command()
@click.argument('targets', nargs=-1)
@click.option('-d', '--dry-run', is_flag=True, default=False)
@click.pass_obj
def run(obj, targets, dry_run):
    """Run the specified workflow."""
    graph = graph_from_config(obj)

    backend_cls = backend_from_config(obj)
    with backend_cls() as backend:
        if targets:
            matched_targets = filter_names(graph, targets)
        else:
            matched_targets = list(graph)

        scheduler = Scheduler(graph=graph, backend=backend, dry_run=dry_run)
        scheduler.schedule_many(matched_targets)
