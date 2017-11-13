import click

from ..backends import backend_from_config
from ..core import graph_from_config
from ..exceptions import TargetNotFoundError


@click.command()
@click.argument('target')
@click.option('-e', '--stderr', is_flag=True)
@click.option('--no-pager', is_flag=True)
@click.pass_obj
def logs(obj, target, stderr, no_pager):
    """Display logs for the latest run of a target.

    By default only standard output is shown. Supply the --stderr flag to show standard error instead.
    """
    graph = graph_from_config(obj)
    backend_cls = backend_from_config(obj)

    if target not in graph:
        raise TargetNotFoundError(target)

    log_file = backend_cls.logs(graph[target], stderr=stderr)
    log_contents = log_file.read()
    log_file.close()

    echo_func = click.echo if no_pager else click.echo_via_pager
    echo_func(log_contents)
