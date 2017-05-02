import statusbar

from ..backends.base import Status
from ..cli import pass_graph, pass_backend
from ..utils import dfs

import click

FILTERS = []


def register_filter(filter_cls):
    FILTERS.append(filter_cls)


class Criteria:
    """A container for filtering criteria."""
    def __init__(self, **kwargs):
        self.__dict__ = kwargs


class FilterType(type):
    def __new__(meta, name, bases, class_dict):
        cls = type.__new__(meta, name, bases, class_dict)
        if cls.__name__ == 'Filter':
            return cls
        register_filter(cls)
        return cls


class Filter(metaclass=FilterType):
    def __init__(self, graph, backend, criteria):
        self.graph = graph
        self.backend = backend
        self.criteria = criteria

    def apply(self, targets):
        return (target for target in targets if self.predicate(target))


class StatusFilter(Filter):
    abstract = True

    def use(self):
        return self.criteria.status

    def predicate(self, target):
        print(target, self.criteria.status, self.backend.status(target))
        if self.criteria.status == 'completed':
            return not self.graph.should_run(target)
        if self.criteria.status == 'shouldrun':
            return self.graph.should_run(target) and self.backend.status(target) == Status.UNKNOWN
        if self.criteria.status == 'running':
            return self.backend.status(target) == Status.RUNNING
        if self.criteria.status == 'submitted':
            return self.backend.status(target) == Status.SUBMITTED


class NameFilter(Filter):
    def use(self):
        return self.criteria.targets

    def predicate(self, target):
        return target.name in self.criteria.targets


class EndpointFilter(Filter):
    def use(self):
        return not self.criteria.all and not self.criteria.targets

    def predicate(self, target):
        return target in self.graph.endpoints()


def filter(graph, backend, criteria):
    targets = iter(graph.targets.values())
    for filter_cls in FILTERS:
        filter = filter_cls(graph, backend, criteria)
        if filter.use():
            targets = filter.apply(targets)
    return targets


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
    filtered_targets = filter(graph, backend, Criteria(**criteria))
    filtered_targets = sorted(filtered_targets, key=lambda t: t.name)

    if names_only:
        for target in filtered_targets:
            click.echo(target.name)
        return

    print_progress(backend, graph, filtered_targets)
