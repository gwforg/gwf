import copy
import logging
import os
import os.path
from pkg_resources import iter_entry_points

import click
from click_plugins import with_plugins

from . import __version__
from .conf import config
from .backends import list_backends
from .exceptions import ConfigurationError

logger = logging.getLogger(__name__)


BASIC_FORMAT = '%(message)s'

ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(message)s'

LOGGING_FORMATS = {
    'warning': BASIC_FORMAT,
    'info': BASIC_FORMAT,
    'debug': ADVANCED_FORMAT,
    'error': BASIC_FORMAT,
}

VERBOSITY_LEVELS = ['warning', 'debug', 'info', 'error']


def ensure_dir(path):
    """Create directory unless it already exists."""
    os.makedirs(path, exist_ok=True)


def get_level(level):
    return getattr(logging, level.upper())


def configure_logging(level_name, formatter_cls):
    fmt = LOGGING_FORMATS[level_name]

    handler = logging.StreamHandler()
    handler.setFormatter(formatter_cls(fmt=fmt))

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(get_level(level_name))


class ColorFormatter(logging.Formatter):

    STYLING = {
        'WARNING': dict(fg='yellow'),
        'INFO': dict(fg='blue'),
        'DEBUG': dict(fg='white'),
        'ERROR': dict(fg='red', bold=True),
        'CRITICAL': dict(fg='magenta', bold=True),
    }

    def format(self, record):
        level = record.levelname
        color_record = copy.copy(record)
        if record.levelname in self.STYLING:
            styling = self.STYLING[level]
            color_record.levelname = click.style(record.levelname, **styling)
            color_record.msg = click.style(record.msg, **styling)
        return super().format(color_record)


@config.validator('backend')
def validate_backend(value):
    backends = list_backends()
    if value not in backends:
        raise ConfigurationError(
            'Invalid value "{}" for key "backend", must be one of: {}.'.format(
                value,
                ', '.join(backends)
            )
        )


@config.validator('verbose')
def validate_verbose(value):
    if value not in VERBOSITY_LEVELS:
        raise ConfigurationError(
            'Invalid value "{}" for key "verbose", must be one of: {}.'.format(
                value,
                ', '.join(VERBOSITY_LEVELS)
            )
        )


@config.validator('no_color')
def validate_no_color(value):
    choices = ('true', 'yes', 'false', 'no')
    if value not in (True, False):
        raise ConfigurationError(
            'Invalid value "{}" for key "no_color", must be one of: {}.'.format(
                value,
                ', '.join(choices)
            )
        )


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
    type=click.Choice(list_backends()),
    default=config['backend'],
    help='Backend used to run workflow.'
)
@click.option(
    '-v',
    '--verbose',
    type=click.Choice(VERBOSITY_LEVELS),
    default=config['verbose'],
    help='Verbosity level.',
)
@click.option(
    '--no-color/--use-color',
    default=config['no_color'],
    help='Enable or disable output colors.'
)
@click.pass_context
def main(ctx, file, backend, verbose, no_color):
    """A flexible, pragmatic workflow tool.

    See help for each command using the `--help` flag for that command:

        gwf status --help

    Shows help for the status command.
    """
    ensure_dir(os.path.join('.gwf'))
    ensure_dir(os.path.join('.gwf', 'logs'))

    formatter_cls = logging.Formatter if no_color else ColorFormatter
    configure_logging(level_name=verbose, formatter_cls=formatter_cls)

    ctx.obj = {'file': file, 'backend': backend}
