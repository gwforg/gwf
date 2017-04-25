import os.path
import logging

import click

from gwf import Graph
from . import __version__
from .utils import WorkflowPath, load_workflow, ensure_dir

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


class WorkflowPathParamType(click.ParamType):
    name = 'workflow'

    def convert(self, value, param, ctx):
        try:
            return WorkflowPath.from_path(value)
        except ValueError:
            self.fail('%s is not a valid integer' % value, param, ctx)

workflow_path = WorkflowPathParamType()


@click.group(context_settings={'obj': {}})
@click.option('-f', '--file', default=WorkflowPath.from_path('workflow.py:gwf'), type=workflow_path, help='Workflow/obj to load')
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
@click.argument('targets', nargs=-1)
@click.pass_context
def run(ctx, targets):
    """Run the specified workflow."""
    graph = Graph(targets=ctx.obj['workflow'].targets)
    print(graph)


@main.command()
@click.argument('targets', nargs=-1)
@click.pass_context
def status(ctx, targets):
    """Show status of a workflow."""
    click.echo('Hi!')


@main.command()
@click.argument('targets', nargs=-1)
@click.pass_context
def clean(ctx, targets):
    """Clean output files of targets."""
    click.echo('Hello')

