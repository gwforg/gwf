import atexit
import logging
import os
import os.path
from functools import update_wrapper

import click
from click_plugins import with_plugins
from pkg_resources import iter_entry_points

from . import __version__
from .config import Config
from .core import Graph
from .utils import ColorFormatter, parse_path, load_workflow, ensure_dir

logger = logging.getLogger(__name__)


BASIC_FORMAT = '%(message)s'

ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(message)s'

LOGGING_FORMATS = {
    'warning': BASIC_FORMAT,
    'info': BASIC_FORMAT,
    'debug': ADVANCED_FORMAT,
    'error': BASIC_FORMAT,
}

BACKENDS = {ep.name: ep.load() for ep in iter_entry_points('gwf.backends')}

CONFIG_DEFAULTS = {
    'verbose': 'info',
    'backend': 'local',
    'no_color': False,
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


def pass_backend(f):
    """Pass the initialized backend to the function."""
    @click.pass_context
    def new_func(ctx, *args, **kwargs):
        backend_name = ctx.obj['_backend']
        backend_cls = BACKENDS[backend_name]

        working_dir = ctx.obj['_workflow'].working_dir
        config = ctx.obj['_config'].get(backend_name, {})

        backend = backend_cls(working_dir=working_dir, config=config)
        atexit.register(backend.close)

        return ctx.invoke(f, *args, backend=backend, **kwargs)
    return update_wrapper(new_func, f)


def pass_graph(f):
    """Pass the complete workflow graph to the function."""
    @click.pass_context
    def new_func(ctx, *args, **kwargs):
        graph = Graph(targets=ctx.obj['_workflow'].targets)
        return ctx.invoke(f, *args, graph=graph, **kwargs)
    return update_wrapper(new_func, f)


def pass_config(f):
    """Pass the complete workflow graph to the function."""
    @click.pass_context
    def new_func(ctx, *args, **kwargs):
        return ctx.invoke(f, *args, config=ctx.obj['_config'], **kwargs)
    return update_wrapper(new_func, f)


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
    help='Backend used to run workflow (default: local).'
)
@click.option(
    '-v',
    '--verbose',
    type=click.Choice(['warning', 'debug', 'info', 'error']),
    help='Verbosity level (default: info).',
)
@click.option(
    '--no-color/--use-color',
    help='Enable or disable output colors.'
)
@click.pass_context
def main(ctx, file, **kwargs):
    """A flexible, pragmatic workflow tool."""
    basedir, filename, obj = parse_path(file)
    workflow = load_workflow(basedir, filename, obj)

    config = Config(
        path=os.path.join(workflow.working_dir, '.gwfconf.json'),
        defaults=CONFIG_DEFAULTS,
        override={k: v for k, v in kwargs.items() if v is not None})

    ensure_dir(os.path.join(workflow.working_dir, '.gwf'))
    ensure_dir(os.path.join(workflow.working_dir, '.gwf', 'logs'))

    backend = config.get('backend')
    verbose = config.get('verbose')
    no_color = config.get('no_color')

    formatter_cls = logging.Formatter if no_color else ColorFormatter
    configure_logging(level_name=verbose, formatter_cls=formatter_cls)

    ctx.obj['_workflow'] = workflow
    ctx.obj['_backend'] = backend
    ctx.obj['_config'] = config
    ctx.obj['working_dir'] = workflow.working_dir
