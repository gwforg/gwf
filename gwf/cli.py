import atexit
import os.path
import logging
from functools import update_wrapper
from pkg_resources import iter_entry_points

import click
from click_plugins import with_plugins

from gwf import Graph
from . import __version__
from .utils import parse_path, load_workflow, ensure_dir


logger = logging.getLogger(__name__)


def get_level(level):
    return getattr(logging, level.upper())

BASIC_FORMAT = '%(message)s'
ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(message)s'

LOGGING_FORMATS = {
    'warning': BASIC_FORMAT,
    'info': BASIC_FORMAT,
    'debug': ADVANCED_FORMAT,
    'error': ADVANCED_FORMAT,
}

BACKENDS = {ep.name: ep.load() for ep in iter_entry_points('gwf.backends')}


def pass_backend(f):
    """Pass the initialized backend to the function."""
    @click.pass_context
    def new_func(ctx, *args, **kwargs):
        backend_cls = BACKENDS[ctx.obj['_backend']]
        backend = backend_cls(working_dir=ctx.obj['_workflow'].working_dir)
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
    default='local',
    type=click.Choice(BACKENDS.keys()),
    help='Backend used to run workflow.'
)
@click.option(
    '-v',
    '--verbose',
    type=click.Choice(['debug', 'info', 'warning', 'error']),
    default='info'
)
@click.pass_context
def main(ctx, backend, file, verbose):
    """A flexible, pragmatic workflow tool."""
    logging.basicConfig(level=get_level(verbose), format=LOGGING_FORMATS[verbose])

    basedir, filename, obj = parse_path(file)
    workflow = load_workflow(basedir, filename, obj)

    ensure_dir(os.path.join(workflow.working_dir, '.gwf'))
    ensure_dir(os.path.join(workflow.working_dir, '.gwf', 'logs'))

    ctx.obj['_workflow'] = workflow
    ctx.obj['_backend'] = backend
