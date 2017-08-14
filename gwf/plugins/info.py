import click

from collections import OrderedDict

from ..core import Target
from ..cli import pass_graph
from ..exceptions import TargetDoesNotExistError


def print_yaml_like(obj, indent=0):
    """Print an object in a YAML-like format."""
    if isinstance(obj, dict):
        if indent > 0:
            print()
        for key, val in obj.items():
            print(' ' * (indent * 2), end='')
            print('{}: '.format(key), end='')
            print_yaml_like(val, indent + 1)
    elif isinstance(obj, list) or isinstance(obj, set):
        if indent > 0:
            print()
        for itm in obj:
            print(' ' * (indent * 2), end='')
            print('- ', end='')
            print_yaml_like(itm, indent=indent + 1)
    elif isinstance(obj, str):
        lines = obj.strip().splitlines()
        if len(lines) == 1:
            print(lines[0])
        else:
            print()
            for line in lines:
                print(' ' * (indent * 2), end='')
                print('|', line)
    elif isinstance(obj, Target):
        print(obj.name)
    else:
        print(obj)


@click.command()
@click.argument('target')
@pass_graph
def info(graph, target):
    """Display information about a target."""

    if target not in graph.targets:
        raise TargetDoesNotExistError(target)
    target = graph.targets[target]

    info_dict = OrderedDict([
        ('name', target.name),
        ('options', target.options),
        ('inputs', target.inputs),
        ('outputs', target.outputs),
        ('spec', target.spec),
        ('dependencies', graph.dependencies[target]),
        ('dependents', graph.dependents[target]),
    ])

    print_yaml_like(info_dict)
