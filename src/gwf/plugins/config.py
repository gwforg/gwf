import click


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
@click.pass_obj
def get(obj, key):
    """Get the value of KEY."""
    click.echo(obj["config"].get(key, "<not set>"))


@config.command()
@click.argument("key")
@click.argument("value")
@click.pass_obj
def set(obj, key, value):
    """Set the value of KEY.

    The key will be created if it does not exist.
    """
    obj["config"][key] = cast_value(value)
    obj["config"].dump()


@config.command()
@click.argument("key")
@click.pass_obj
def unset(obj, key):
    """Unset KEY."""
    del obj["config"][key]
    obj["config"].dump()
