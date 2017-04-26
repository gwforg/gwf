import statusbar

from ..cli import pass_graph, pass_backend
from ..utils import dfs

import click


def _split_target_list(backend, graph, targets):
    should_run, submitted, running, completed = [], [], [], []
    for target in targets:
        if graph.should_run(target):
            if backend.running(target):
                running.append(target)
            elif backend.submitted(target):
                submitted.append(target)
            else:
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


@click.command()
@click.argument('targets', nargs=-1)
@pass_graph
@pass_backend
def status(backend, graph, targets):
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
    matched_targets = graph.get_targets_by_name(targets) or graph.endpoints()
    sorted_targets = sorted(matched_targets, key=lambda t: t.name)
    print_progress(backend, graph, sorted_targets)
