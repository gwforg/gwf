import statusbar

from ..backends.base import Status
from ..cli import pass_graph, pass_backend
from ..filtering import Criteria, filter
from ..utils import dfs

import click


def _split_target_list(backend, graph, targets):
    should_run, submitted, running, completed = [], [], [], []
    for target in targets:
        status = backend.status(target)
        if status == Status.RUNNING:
            running.append(target)
        elif status == Status.SUBMITTED:
            submitted.append(target)
        elif status == Status.UNKNOWN:
            if graph.should_run(target):
                should_run.append(target)
            else:
                completed.append(target)
    return should_run, submitted, running, completed


def print_progress(backend, graph, targets):
    table = statusbar.StatusTable(fill_char=' ')
    for target in targets:
        dependencies = dfs(target, graph.dependencies)
        should_run, submitted, running, completed = _split_target_list(backend, graph, dependencies)
        status_bar = table.add_status_line(target.name)
        status_bar.add_progress(len(completed), 'C', color='green')
        status_bar.add_progress(len(running), 'R', color='blue')
        status_bar.add_progress(len(submitted), 'S', color='yellow')
        status_bar.add_progress(len(should_run), '.', color='magenta')
    print('\n'.join(table.format_table()))


def _status(backend, graph, names_only, **criteria):
    filtered_targets = filter(graph, backend, Criteria(**criteria))
    filtered_targets = sorted(filtered_targets, key=lambda t: t.name)

    if names_only:
        for target in filtered_targets:
            click.echo(target.name)
        return

    print_progress(backend, graph, filtered_targets)


@click.command()
@click.argument('targets', nargs=-1)
@click.option('-n', '--names-only', is_flag=True)
@click.option('--all/--endpoints', help='Whether to show all targets or only endpoints if no targets are specified.')
@click.option('-s', '--status', type=click.Choice(['shouldrun', 'submitted', 'running', 'completed']))
@pass_graph
@pass_backend
def status(backend, graph, names_only, **criteria):
    """
    Show the status of targets.

    By default, shows a progress bar for each endpoint in the workflow.
    If one or more target names are supplied, progress bars are shown
    for these targets.

    A progress bar represents the target and its dependencies, and
    shows how many of the dependencies either should run (magenta, .),
    are submitted (yellow, S), are running (blue, R), are
    completed (green, C), or have failed (red, F).
    """
    _status(backend, graph, names_only, **criteria)
