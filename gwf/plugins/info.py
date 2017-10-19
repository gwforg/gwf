import json

import click

from collections import OrderedDict

from ..core import graph_from_config
from ..filtering import filter_names


@click.command()
@click.argument('targets', nargs=-1)
@click.pass_obj
def info(obj, targets):
    """Display information about a target."""
    graph = graph_from_config(obj)

    matches = iter(graph)
    if targets:
        matches = filter_names(matches, targets)

    obj = {}
    for target in matches:
        obj[target.name] = OrderedDict([
            ('options', target.options),
            ('inputs', target.inputs),
            ('outputs', target.outputs),
            ('spec', target.spec),
            ('dependencies', [target.name for target in graph.dependencies[target]]),
            ('dependents', [target.name for target in graph.dependents[target]]),
        ])

    print(json.dumps(obj, indent=4))
