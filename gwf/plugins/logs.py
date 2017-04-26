import click

from ..cli import pass_graph, pass_backend
from ..exceptions import TargetDoesNotExistError


@click.command()
@click.argument('target')
@click.option('-o/-e', '--stdout/--stderr')
@pass_graph
@pass_backend
def logs(backend, graph, target, stdout):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show standard error instead.
    """
    if target not in graph.targets:
        raise TargetDoesNotExistError(target)

    log = backend.logs(graph.targets[target], stderr=not stdout)
    click.echo(log.read())
