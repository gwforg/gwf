import json
from collections import OrderedDict

import click

from ..core import Graph
from ..filtering import filter_names
from ..workflow import Workflow


def print_json(targets, graph):
    obj = {}
    for target in targets:
        obj[target.name] = OrderedDict(
            [
                ("options", target.options),
                ("inputs", target.inputs),
                ("outputs", target.outputs),
                ("spec", target.spec),
                (
                    "dependencies",
                    [target.name for target in graph.dependencies[target]],
                ),
                ("dependents", [target.name for target in graph.dependents[target]]),
            ]
        )

    print(json.dumps(obj, indent=4))


FORMATS = {
    "json": print_json,
}


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-f", "--format", type=click.Choice(["json"]), default="json")
@click.pass_obj
def info(obj, targets, format):
    """Display information about a target."""
    workflow = Workflow.from_config(obj)
    graph = Graph.from_targets(workflow.targets)

    matches = iter(graph)
    if targets:
        matches = filter_names(matches, targets)

    printer = FORMATS[format]
    printer(matches, graph)
