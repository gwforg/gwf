import click

from ..core import pass_context


@click.group()
def config():
    """Set, unset and get configuration."""


@config.command()
@click.argument("key")
@pass_context
def get(ctx, key):
    """Get the value of KEY."""
    click.echo(ctx.config.get(key, "<not set>"))


@config.command()
@click.argument("key")
@click.argument("value")
@pass_context
def set(ctx, key, value):
    """Set the value of KEY.

    The key will be created if it does not exist.
    """
    ctx.config[key] = value
    ctx.config.dump()


@config.command()
@click.argument("key")
@pass_context
def unset(ctx, key):
    """Unset KEY."""
    del ctx.config[key]
    ctx.config.dump()
