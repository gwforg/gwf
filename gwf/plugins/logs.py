import click

from ..cli import pass_graph, pass_backend
from ..exceptions import TargetDoesNotExistError


@click.command()
@click.argument('target')
@click.option('-e', '--stderr', is_flag=True)
@click.option('--no-pager', is_flag=True)
@pass_graph
@pass_backend
def logs(backend, graph, target, stderr, no_pager):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show standard error instead.
    """
    if target not in graph.targets:
        raise TargetDoesNotExistError(target)

    log_file = backend.logs(graph.targets[target], stderr=stderr)
    log_contents = log_file.read()
    log_file.close()

    echo_func = click.echo if no_pager else click.echo_via_pager
    echo_func(log_contents)
