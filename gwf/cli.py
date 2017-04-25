import os.path
import logging

import click

from gwf import Graph
from . import __version__
from .exceptions import GWFError
from .utils import load_workflow, ensure_dir

logger = logging.getLogger(__name__)

VERBOSITY_LEVELS = {
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG,
}

BASIC_FORMAT = '%(message)s'
ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(message)s'

LOGGING_FORMATS = {
    'warning': BASIC_FORMAT,
    'info': BASIC_FORMAT,
    'debug': ADVANCED_FORMAT,
    'error': ADVANCED_FORMAT,
}


class WorkflowPath:
    def __init__(self, basedir, filename, obj='gwf'):
        self.basedir = basedir
        self.filename = filename
        self.obj = obj

    @classmethod
    def from_path(cls, path):
        comps = path.rsplit(':')
        if len(comps) == 2:
            path, obj = comps
        elif len(comps) == 1:
            path, obj = comps[0], None
        else:
            raise ValueError('Invalid path: "{}".'.format(path))

        basedir, filename = os.path.split(path)
        filename, _ = os.path.splitext(filename)
        if obj is None:
            return cls(basedir, filename)
        return cls(basedir, filename, obj)


class WorkflowPathParamType(click.ParamType):
    name = 'workflow'

    def convert(self, value, param, ctx):
        try:
            return WorkflowPath.from_path(value)
        except ValueError:
            self.fail('%s is not a valid integer' % value, param, ctx)

workflow_path = WorkflowPathParamType()


@click.group(context_settings={'obj': {}})
@click.option('-f', '--file', type=workflow_path, default=WorkflowPath.from_path('workflow.py:gwf'), help='Workflow/obj to load')
@click.option('-b', '--backend', default='local', help='Backend used to run workflow')
@click.option('-v', '--verbose', type=click.Choice(['debug', 'info', 'warning', 'error']), default='info')
@click.version_option(version=__version__)
@click.pass_context
def main(ctx, backend, file, verbose):
    """A flexible, pragmatic workflow tool."""
    logging.basicConfig(level=VERBOSITY_LEVELS[verbose], format=LOGGING_FORMATS[verbose])

    workflow = load_workflow(file.basedir, file.filename, file.obj)

    ensure_dir(os.path.join(workflow.working_dir, '.gwf'))
    ensure_dir(os.path.join(workflow.working_dir, '.gwf', 'logs'))

    ctx.obj['workflow'] = workflow
    ctx.obj['backend'] = backend


@main.command()
@click.pass_context
def run(ctx):
    """Run the specified workflow."""
    graph = Graph(targets=ctx.obj['workflow'].targets)
    print(graph)


@main.command()
@click.pass_context
def status():
    """Show status of a workflow."""
    click.echo('Hi!')


@main.command()
@click.pass_context
def clean():
    """Clean output files of targets."""
    click.echo('Hello')

