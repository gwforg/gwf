import click

from ..conf import file_config


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
        except Exception:
            continue
    return value


@click.group()
def config():
    """Set, unset and get configuration."""


@config.command()
@click.argument('key')
def get(key):
    """Get the value of KEY."""
    click.echo(file_config.get(key))


@config.command()
@click.argument('key')
@click.argument('value')
def set(key, value):
    """Set the value of KEY.

    The key will be created if it does not exist.
    """
    file_config[key] = cast_value(value)
    file_config.dump()


@config.command()
@click.argument('key')
def unset(key):
    """Unset KEY."""
    del file_config[key]
    file_config.dump()
