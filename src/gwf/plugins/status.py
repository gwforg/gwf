import logging
import os.path
from collections import Counter

import attr
import click

from gwf.backends.base import Status

from ..backends import Backend
from ..conf import config
from ..core import CachedFilesystem, Graph, TargetStatus
from ..filtering import EndpointFilter, NameFilter, StatusFilter, filter_generic
from ..scheduling import schedule_workflow
from ..utils import PersistableDict
from ..workflow import Workflow

logger = logging.getLogger(__name__)


def _status_str_to_enum(s):
    return TargetStatus[s.upper()]


def _status_strs_to_enums(iterable):
    return list(map(_status_str_to_enum, iterable))


@attr.define(frozen=True)
class HumanStatus:
    status: TargetStatus = attr.ib()
    reason: str = attr.ib(eq=False)
    symbol: str = attr.ib()
    color: str = attr.ib()


def get_status(reason, backend):
    """Return the status of a target.

    Returns the status of a target where it is taken into account whether
    the target should run or not.

    :param Target target:
        The target to return status for.
    """
    status = backend.status(reason.target)
    if status == Status.RUNNING:
        return HumanStatus(TargetStatus.RUNNING, "is running", symbol="↻", color="blue")
    elif status == Status.SUBMITTED:
        return HumanStatus(
            TargetStatus.SUBMITTED, "has been submitted", symbol="-", color="cyan"
        )
    elif reason.scheduled:
        return HumanStatus(TargetStatus.SHOULDRUN, reason, symbol="⨯", color="magenta")
    else:
        return HumanStatus(TargetStatus.COMPLETED, reason, symbol="✓", color="green")


def print_table(graph, targets, target_status):
    targets = list(targets)

    name_col_width = max((len(target.name) for target in targets), default=0) + 4
    format_str = "{status} {name:<{name_col_width}}{percentage:>7.2%}    {explanation}"

    for target in sorted(targets, key=lambda t: t.order):
        status = target_status[target]
        deps = graph.dfs(target)
        percentage = sum(
            1 for t in deps if target_status[t].status == TargetStatus.COMPLETED
        ) / len(deps)
        line = format_str.format(
            name=target.name,
            status=status.symbol,
            percentage=percentage,
            name_col_width=name_col_width,
            explanation=status.reason or "",
        )
        click.secho(line, fg=status.color)


def print_summary(_, targets, target_status):
    status_counts = Counter(status for status in target_status.values())
    click.echo("∑ {:<10}{:>10}".format("total", len(targets)))
    for status, count in status_counts.items():
        click.secho(
            "{} {:<10}{:>10}".format(
                status.symbol,
                status.status.name.lower(),
                count,
            ),
            fg=status.color,
        )


FORMATS = {
    "summary": print_summary,
    "default": print_table,
}


@click.command()
@click.argument("targets", nargs=-1)
@click.option("--endpoints", is_flag=True, default=False, help="Show only endpoints.")
@click.option(
    "-f",
    "--format",
    default="default",
    type=click.Choice(["summary", "default"]),
    help="How to format status output.",
)
@click.option(
    "-s",
    "--status",
    type=click.Choice(["shouldrun", "submitted", "running", "completed"]),
    multiple=True,
)
@click.pass_obj
def status(obj, status, endpoints, format, targets):
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

    spec_hashes = None
    if config.get("use_spec_hashes"):
        spec_hashes = PersistableDict(
            os.path.join(workflow.working_dir, ".gwf", "spec-hashes.json")
        )

    with backend_cls() as backend:
        reasons = schedule_workflow(
            graph, fs=CachedFilesystem(), spec_hashes=spec_hashes
        )

        target_status = {}
        for target, reason in reasons.items():
            target_status[target] = get_status(reason, backend)

        def status_provider(target):
            return target_status[target].status

        filters = []
        if status:
            filters.append(
                StatusFilter(
                    status_provider=status_provider,
                    status=_status_strs_to_enums(status),
                )
            )
        if targets:
            filters.append(NameFilter(patterns=targets))
        if endpoints:
            filters.append(EndpointFilter(endpoints=graph.endpoints()))

        matches = filter_generic(targets=graph, filters=filters)

        printer = FORMATS[format]
        printer(graph, matches, target_status)
