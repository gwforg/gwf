import click

from ..conf import config as _config


def humanbool(x):
    if x in ("true", "yes"):
        return True
    elif x in ("false", "no"):
        return False
    raise TypeError("x is not a boolean.")


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
@click.argument("key")
def get(key):
    """Get the value of KEY."""
    click.echo(_config.get(key, "<not set>"))


@config.command()
@click.argument("key")
@click.argument("value")
def set(key, value):
    """Set the value of KEY.

    The key will be created if it does not exist.
    """
    _config[key] = cast_value(value)
    _config.dump()


@config.command()
@click.argument("key")
def unset(key):
    """Unset KEY."""
    del _config[key]
    _config.dump()
