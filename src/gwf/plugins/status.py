import click

from ..backends import Status, backend_from_config
from ..core import Scheduler, graph_from_config
from ..filtering import StatusFilter, EndpointFilter, NameFilter, filter_generic


STATUS_COLORS = {
    'UNKNOWN': 'red',
    'SHOULD RUN': 'red',
    'RUNNING': 'blue',
    'COMPLETED': 'green',
}




def print_table(backend, graph, targets):
    scheduler = Scheduler(backend=backend, graph=graph)

    name_col_width = max((len(target.name) for target in targets), default=0) + 1
    format_str = '{name:<{name_col_width}}\t{status:<10}'

    for target in targets:
        status = backend.status(target)
        if status == Status.UNKNOWN:
            if scheduler.should_run(target):
                status = 'SHOULD RUN'
            else:
                status = 'COMPLETED'
        else:
            status = status.name

        color = STATUS_COLORS[status]

        line = format_str.format(
            name=target.name,
            status=click.style(status, fg=color),
            name_col_width=name_col_width
        )
        click.echo(line)


@click.command()
@click.argument(
    'targets', nargs=-1
)
@click.option(
    '--all',
    is_flag=True,
    default=False,
    help='Show all targets, not only endpoints.'
)
@click.option(
    '-s', '--status',
    type=click.Choice(['shouldrun', 'submitted', 'running', 'completed'])
)
@click.pass_obj
def status(obj, status, all, targets):
    """
    Show the status of targets.

    By default, shows a progress bar for each endpoint in the workflow.
    If one or more target names are supplied, progress bars are shown
    for these targets.

    A progress bar represents the target and its dependencies, and
    shows how many of the dependencies either should run (magenta, .),
    are submitted (yellow, S), are running (blue, R), or are
    completed (green, C).
    """
    graph = graph_from_config(obj)
    backend_cls = backend_from_config(obj)

    with backend_cls() as backend:
        scheduler = Scheduler(graph=graph, backend=backend)

        filters = []
        if status:
            filters.append(StatusFilter(scheduler=scheduler, status=status))
        if targets:
            filters.append(NameFilter(patterns=targets))
        if not all:
            filters.append(EndpointFilter(endpoints=graph.endpoints()))

        matches = filter_generic(targets=graph, filters=filters)
        matches = sorted(matches, key=lambda t: t.name)
        print_table(backend, graph, matches)
