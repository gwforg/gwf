from ..core import schedule_many
from ..cli import pass_backend, pass_graph

import click


@click.command()
@click.argument('targets', nargs=-1)
@click.option('-d', '--dry-run', is_flag=True, default=False)
@pass_backend
@pass_graph
def run(graph, backend, targets, dry_run):
    """Run the specified workflow."""
    matched_targets = graph.get_targets_by_name(targets) or graph.endpoints()
    schedule_many(graph, backend, matched_targets, dry_run=dry_run)
