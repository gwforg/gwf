import logging
import os
import os.path
from pkg_resources import iter_entry_points

import click
from click_plugins import with_plugins

from . import __version__
from .conf import config
from .backends import BACKENDS
from .utils import ColorFormatter, ensure_dir

logger = logging.getLogger(__name__)


BASIC_FORMAT = '%(message)s'

ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(message)s'

LOGGING_FORMATS = {
    'warning': BASIC_FORMAT,
    'info': BASIC_FORMAT,
    'debug': ADVANCED_FORMAT,
    'error': BASIC_FORMAT,
}


def get_level(level):
    return getattr(logging, level.upper())


def configure_logging(level_name, formatter_cls):
    fmt = LOGGING_FORMATS[level_name]

    handler = logging.StreamHandler()
    handler.setFormatter(formatter_cls(fmt=fmt))

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(get_level(level_name))


@with_plugins(iter_entry_points('gwf.plugins'))
@click.group(context_settings={'obj': {}})
@click.version_option(version=__version__)
@click.option(
    '-f',
    '--file',
    default='workflow.py:gwf',
    help='Workflow/obj to load.'
)
@click.option(
    '-b',
    '--backend',
    type=click.Choice(BACKENDS.keys()),
    default=config.get('backend', 'local'),
    help='Backend used to run workflow.'
)
@click.option(
    '-v',
    '--verbose',
    type=click.Choice(['warning', 'debug', 'info', 'error']),
    default=config.get('verbose', 'info'),
    help='Verbosity level.',
)
@click.option(
    '--no-color/--use-color',
    default=config.get('no_color', False),
    help='Enable or disable output colors.'
)
@click.pass_context
def main(ctx, file, backend, verbose, no_color):
    """A flexible, pragmatic workflow tool."""
    ensure_dir(os.path.join('.gwf'))
    ensure_dir(os.path.join('.gwf', 'logs'))

    formatter_cls = logging.Formatter if no_color else ColorFormatter
    configure_logging(level_name=verbose, formatter_cls=formatter_cls)

    ctx.obj = {'file': file, 'backend': backend}
