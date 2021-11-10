from collections import Counter

import click

from ..backends import Backend
from ..core import Graph, TargetStatus, get_status, schedule
from ..filtering import EndpointFilter, NameFilter, StatusFilter, filter_generic
from ..workflow import Workflow

STATUS_COLORS = {
    TargetStatus.SHOULDRUN: "magenta",
    TargetStatus.SUBMITTED: "yellow",
    TargetStatus.RUNNING: "blue",
    TargetStatus.COMPLETED: "green",
}

STATUS_ORDER = (
    TargetStatus.SHOULDRUN,
    TargetStatus.SUBMITTED,
    TargetStatus.RUNNING,
    TargetStatus.COMPLETED,
)


def _status_str_to_enum(s):
    return TargetStatus[s.upper()]


def _status_strs_to_enums(iterable):
    return list(map(_status_str_to_enum, iterable))


def print_table(graph, targets, status_provider):
    targets = list(targets)

    name_col_width = max((len(target.name) for target in targets), default=0) + 4
    format_str = (
        "{name:<{name_col_width}}{status:<23}{percentage:>7.2%}"
        " [{num_shouldrun}/{num_submitted}/{num_running}/{num_completed}]"
    )

    for target in sorted(targets, key=lambda t: t.order):
        status = status_provider(target)
        color = STATUS_COLORS[status]

        deps = graph.dfs(target)
        deps_total = len(deps)

        def num_deps_with_status(status):
            return sum(1 for target in deps if status_provider(target) == status)

        num_shouldrun = num_deps_with_status(TargetStatus.SHOULDRUN)
        num_submitted = num_deps_with_status(TargetStatus.SUBMITTED)
        num_running = num_deps_with_status(TargetStatus.RUNNING)
        num_completed = num_deps_with_status(TargetStatus.COMPLETED)

        percentage = num_completed / deps_total

        line = format_str.format(
            name=target.name,
            status=click.style(status.name.lower(), fg=color),
            percentage=percentage,
            num_shouldrun=num_shouldrun,
            num_submitted=num_submitted,
            num_running=num_running,
            num_completed=num_completed,
            name_col_width=name_col_width,
        )
        click.echo(line)


def print_summary(backend, graph, targets):
    status_counts = Counter(get_status(target, backend=backend) for target in targets)
    click.echo("{:<15}{:>10}".format("total", len(targets)))
    for status in STATUS_ORDER:
        color = STATUS_COLORS[status]
        padded_name = "{:<15}".format(status.name.lower())
        click.echo(
            "{}{:>10}".format(click.style(padded_name, fg=color), status_counts[status])
        )


@click.command()
@click.argument("targets", nargs=-1)
@click.option("--endpoints", is_flag=True, default=False, help="Show only endpoints.")
@click.option(
    "--summary", is_flag=True, default=False, help="Only show summary statistics."
)
@click.option(
    "-s",
    "--status",
    type=click.Choice(["shouldrun", "submitted", "running", "completed"]),
    multiple=True,
)
@click.pass_obj
def status(obj, status, summary, endpoints, targets):
    """
    Show the status of targets.

    One target per line is shown. Each line contains the target name, the
    status of the target itself, and the percentage of dependencies of the
    target that are completed (including the target itself). That is, 100%
    means that the target and all of its dependencies have been completed.

    In square brackets, the number of targets that should run, have been
    submitted, are running, and completed targets are shown, respectively.

    The `-s/--status` flag can be applied multiple times to show targets that
    match either of the queries, e.g. `gwf status -s shouldrun -s completed`
    will show all targets that should run and all targets that are completed.

    The targets are shown in creation-order.
    """
    workflow = Workflow.from_config(obj)
    graph = Graph.from_targets(workflow.targets)
    backend_cls = Backend.from_config(obj)

    scheduled, _ = schedule(graph.endpoints(), graph=graph)

    def status_provider(target):
        return get_status(target, scheduled, backend)

    with backend_cls() as backend:
        filters = []
        if status:
            status = _status_strs_to_enums(status)
            filters.append(StatusFilter(status_provider=status_provider, status=status))
        if targets:
            filters.append(NameFilter(patterns=targets))
        if endpoints:
            filters.append(EndpointFilter(endpoints=graph.endpoints()))

        matches = filter_generic(targets=graph, filters=filters)

        if not summary:
            print_table(graph, matches, status_provider)
        else:
            print_summary(graph, matches, status_provider)
