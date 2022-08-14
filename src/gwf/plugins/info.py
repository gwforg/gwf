import json
from collections import OrderedDict

import click

from ..core import Graph
from ..filtering import filter_names
from ..workflow import Workflow


def _indent(s):
    return "    " + s


def print_label(label):
    click.secho(label, bold=True)


def print_value(value, as_filename=False):
    if as_filename:
        value = click.format_filename(value)
    click.secho(_indent(value))


def print_list(values, **kwargs):
    if not values:
        print_value("-")
    for value in values:
        print_value(value, **kwargs)


def print_pretty(targets, graph):
    for target in targets:
        print_label("Name:")
        print_value(target.name)
        print_label("Inputs:")
        print_list(target.inputs, as_filename=True)
        print_label("Outputs:")
        print_list(target.outputs, as_filename=True)
        print_label("Dependents:")
        print_list([target.name for target in graph.dependents[target]])
        print_label("Spec:")
        print_list(target.spec.strip().splitlines())
        click.secho()


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
    "pretty": print_pretty,
}


@click.command()
@click.argument("targets", nargs=-1)
@click.option("-f", "--format", type=click.Choice(["json", "pretty"]), default="json")
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
