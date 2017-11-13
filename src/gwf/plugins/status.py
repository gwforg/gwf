import click
import statusbar

from ..backends import Status, backend_from_config
from ..core import Scheduler, graph_from_config
from ..filtering import StatusFilter, EndpointFilter, NameFilter, filter_generic


def dfs(root, dependencies):
    """Return the depth-first traversal path through a graph from `root`."""
    visited = set()
    path = []

    def dfs_inner(node):
        if node in visited:
            return

        visited.add(node)
        for dep in dependencies[node]:
            dfs_inner(dep)

        path.append(node)

    dfs_inner(root)
    return path


def split_target_list(backend, scheduler, targets):
    should_run, submitted, running, completed = [], [], [], []
    for target in targets:
        status = backend.status(target)
        if status == Status.RUNNING:
            running.append(target)
        elif status == Status.SUBMITTED:
            submitted.append(target)
        elif status == Status.UNKNOWN:
            if scheduler.should_run(target):
                should_run.append(target)
            else:
                completed.append(target)
    return should_run, submitted, running, completed


def print_progress(backend, graph, targets):
    scheduler = Scheduler(graph=graph, backend=backend)
    table = statusbar.StatusTable(fill_char=' ')
    for target in targets:
        dependencies = dfs(target, graph.dependencies)
        should_run, submitted, running, completed = split_target_list(backend, scheduler, dependencies)
        status_bar = table.add_status_line(target.name)
        status_bar.add_progress(len(completed), 'C', color='green')
        status_bar.add_progress(len(running), 'R', color='blue')
        status_bar.add_progress(len(submitted), 'S', color='yellow')
        status_bar.add_progress(len(should_run), '.', color='magenta')
    click.echo('\n'.join(table.format_table()))


def print_table(backend, graph, targets):
    scheduler = Scheduler(backend=backend, graph=graph)
    
    name_col_width = max((len(target.name) for target in targets), default=0) + 1
    format_str = '{name:<{name_col_width}}{status:<10}'

    for target in targets:
        status = backend.status(target)
        if status == Status.UNKNOWN:
            if scheduler.should_run(target):
                status = 'SHOULDRUN'
            else:
                status = 'COMPLETED'
        else:
            status = status.name

        click.echo(format_str.format(name=target.name, status=status, name_col_width=name_col_width))


@click.command()
@click.argument(
    'targets', nargs=-1
)
@click.option(
    '--format',
    type=click.Choice(['progress', 'table']),
    default='progress'
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
def status(obj, status, all, format, targets):
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
    format_funcs = {
        'progress': print_progress,
        'table': print_table,
    }

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
        format_func = format_funcs[format]
        format_func(backend, graph, matches)
