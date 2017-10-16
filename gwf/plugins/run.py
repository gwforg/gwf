from ..core import Scheduler
from ..cli import pass_backend, pass_graph
from ..filtering import filter_names

import click


@click.command()
@click.argument('targets', nargs=-1)
@click.option('-d', '--dry-run', is_flag=True, default=False)
@pass_backend
@pass_graph
def run(graph, backend, targets, dry_run):
    """Run the specified workflow."""
    matched_targets = filter_names(graph.targets.values(), targets)
    scheduler = Scheduler(graph=graph, backend=backend, dry_run=dry_run)
    scheduler.schedule_many(matched_targets)
