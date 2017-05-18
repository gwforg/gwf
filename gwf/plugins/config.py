import json

import click

from ..cli import pass_config


def humanbool(x):
    if x in ('true', 'yes'):
        return True
    elif x in ('false', 'no'):
        return False
    raise TypeError('x is not a boolean.')


def cast_value(value):
    types = [int, humanbool, str]
    for type_func in types:
        try:
            return type_func(value)
        except:
            continue
    return value


@click.group()
def config():
    """Set, unset and retrieve configuration."""


@config.command()
@click.argument('key')
@pass_config
def get(config, key):
    """Get the value of KEY."""
    value = config.get(key)
    click.echo(json.dumps(value, indent=4, sort_keys=True))


@config.command()
@click.argument('key')
@click.argument('value')
@pass_config
def set(config, key, value):
    """Set the value of KEY.

    The key will be created if it does not exist.
    """
    file_config = config.get_file_config()
    file_config.set(key, cast_value(value))
    file_config.dump()
