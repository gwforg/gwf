import logging
from collections import Counter, defaultdict

import click

from .. import Workflow
from ..backends import create_backend
from ..core import CachedFilesystem, Graph, Status, get_spec_hashes, pass_context
from ..filtering import (
    EndpointFilter,
    GroupFilter,
    NameFilter,
    StatusFilter,
    filter_generic,
)
from ..scheduling import get_status_map

logger = logging.getLogger(__name__)


_STATUS_VISUALS = {
    Status.SHOULDRUN: ("magenta", "⨯"),
    Status.SUBMITTED: ("cyan", "-"),
    Status.RUNNING: ("blue", "↻"),
    Status.COMPLETED: ("green", "✓"),
    Status.FAILED: ("red", "⨯"),
    Status.CANCELLED: ("red", "⨯"),
}

_STATUS_NAMES = [s.name.lower() for s in Status]


def _status_names_to_enums(iterable):
    return [Status[s.upper()] for s in iterable]


def _key_target_order(item):
    t, _ = item
    return t.order


def print_grouped_table(target_states, backend):
    grouped = defaultdict(list)
    for target, target_state in target_states.items():
        grouped[target.group or "-"].append((target, target_state))

    name_col_width = max((len(group) for group in grouped), default=0) + 4
    for group, target_states in sorted(grouped.items()):
        click.secho(group.ljust(name_col_width), nl=False)
        status_counts = Counter(status for _, status in target_states)
        _, max_count = status_counts.most_common()[0]
        count_width = len(str(max_count)) + 2
        for status, (color, _) in _STATUS_VISUALS.items():
            click.secho(
                f"{status_counts[status]:>{count_width}} {status.name.lower():<10}",
                fg=color,
                nl=False,
            )
        click.secho()


def print_flat_table(target_states, backend):
    name_col_width = (
        max((len(target.name) for target in target_states.values()), default=0) + 4
    )
    format_str = "{symbol} {name:<{name_col_width}} {status} (id: {tracked_id})"
    for target, status in sorted(target_states.items(), key=_key_target_order):
        color, symbol = _STATUS_VISUALS[status]

        if hasattr(backend, "get_tracked_id"):
            tracked_id = backend.get_tracked_id(target) or "none"
        else:
            tracked_id = "none"

        line = format_str.format(
            symbol=symbol,
            name=target.name,
            status=status.name.lower(),
            name_col_width=name_col_width,
            tracked_id=tracked_id,
        )
        click.secho(line, fg=color)


def print_summary(target_states, backend):
    status_counts = Counter(status for status in target_states.values())
    _, max_count = status_counts.most_common()[0]
    count_width = len(str(max_count)) + 4
    for status, (color, symbol) in _STATUS_VISUALS.items():
        click.secho(
            f"{symbol} {status.name.lower():<10}{status_counts[status]:>{count_width}}",
            fg=color,
        )


FORMATS = {
    "summary": print_summary,
    "default": print_flat_table,
    "grouped": print_grouped_table,
}


@click.command()
@click.argument("targets", nargs=-1)
@click.option("--endpoints", is_flag=True, default=False, help="Show only endpoints.")
@click.option(
    "-f",
    "--format",
    default="default",
    type=click.Choice(FORMATS.keys()),
    help="How to format status output.",
)
@click.option(
    "-s",
    "--status",
    type=click.Choice(_STATUS_NAMES),
    multiple=True,
)
@click.option(
    "-g",
    "--group",
    multiple=True,
)
@pass_context
def status(ctx, group, status, endpoints, format, targets):
    """
    Show the status of targets.

    With the default format, one target per line is shown. Each line contains
    the target name, the status of the target itself, and the id that identifies
    the target in the backend (for example, the job id in the case of Slurm).

    The grouped format will show targets grouped according to their `group`
    parameter. Targets where `group=None` will be grouped into the `-` group.
    For each group the number of targets in each state is shown.

    The summary format will show a short summary of the target states for the
    entire workflow.

    The `-s/--status` flag can be applied multiple times to show targets that
    match either of the queries, e.g. `gwf status -s shouldrun -s completed`
    will show all targets that should run and all targets that are completed.
    When the grouped format is used, groups containing at least one target
    matching the filter will be shown.

    The targets are shown in creation-order.
    """
    workflow = Workflow.from_context(ctx)

    fs = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, fs)

    with create_backend(
        ctx.backend, working_dir=ctx.working_dir, config=ctx.config
    ) as backend:
        with get_spec_hashes(
            working_dir=ctx.working_dir, config=ctx.config
        ) as spec_hashes:
            target_states = get_status_map(
                graph,
                fs,
                spec_hashes,
                backend,
            )

            filters = []
            if status:
                filters.append(
                    StatusFilter(
                        status_provider=target_states.get,
                        status=_status_names_to_enums(status),
                    )
                )
            if targets:
                filters.append(NameFilter(patterns=targets))
            if endpoints:
                filters.append(EndpointFilter(endpoints=graph.endpoints()))
            if group:
                filters.append(GroupFilter(patterns=group))

            matches = set(filter_generic(targets=graph, filters=filters))
            target_states = {k: v for k, v in target_states.items() if k in matches}

            printer = FORMATS[format]
            printer(target_states, backend)
