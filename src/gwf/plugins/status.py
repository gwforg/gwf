import click

from ..backends import Status, backend_from_config
from ..core import Scheduler, graph_from_config
from ..filtering import StatusFilter, EndpointFilter, NameFilter, filter_generic


STATUS_COLORS = {
    'submitted': 'yellow',
    'shouldrun': 'magenta',
    'running': 'blue',
    'completed': 'green',
}


def _get_status(target, backend, scheduler):
    status = backend.status(target)
    if status == Status.UNKNOWN:
        if scheduler.should_run(target):
            return 'shouldrun'
        else:
            return 'completed'
    return status.name.lower()


def print_table(backend, graph, targets):
    scheduler = Scheduler(backend=backend, graph=graph)

    name_col_width = max((len(target.name) for target in targets), default=0) + 4
    format_str = '{name:<{name_col_width}}{status:<23}{percentage:>7.2%}'

    for target in targets:
        status = _get_status(target, backend, scheduler)
        color = STATUS_COLORS[status]

        deps = graph.dfs(target)
        deps_total = len(deps)
        deps_completed = len([
            target
            for target in deps 
            if _get_status(target, backend, scheduler) == 'completed'
        ])
        percentage = deps_completed / deps_total

        line = format_str.format(
            name=target.name,
            status=click.style(status, fg=color),
            percentage=percentage,
            name_col_width=name_col_width
        )
        click.echo(line)


@click.command()
@click.argument(
    'targets', nargs=-1
)
@click.option(
    '--endpoints',
    is_flag=True,
    default=False,
    help='Show only endpoints.'
)
@click.option(
    '-s', '--status',
    type=click.Choice(['shouldrun', 'submitted', 'running', 'completed'])
)
@click.pass_obj
def status(obj, status, endpoints, targets):
    """
    Show the status of targets.

    One target per line is shown. Each line contains the target name, the status
    of the target itself, and the percentage of dependencies of the target that
    are completed (including the target itself). That is, 100% means that the
    target and all of its dependencies have been completed.

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
        if endpoints:
            filters.append(EndpointFilter(endpoints=graph.endpoints()))

        matches = filter_generic(targets=graph, filters=filters)
        matches = sorted(matches, key=lambda t: t.name)
        print_table(backend, graph, matches)
