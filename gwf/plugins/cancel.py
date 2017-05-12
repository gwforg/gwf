from ..exceptions import UnsupportedOperationError
from ..backends.base import UnknownTargetError
from ..cli import pass_backend, pass_graph

import click


@click.command()
@click.argument('targets', nargs=-1)
@pass_backend
@pass_graph
def cancel(graph, backend, targets):
    """Cancel the specified targets."""
    matched_targets = graph.find(targets)
    if not matched_targets:
        if click.confirm('This will cancel all targets! Do you want to continue?', abort=True):
            matched_targets = graph.targets.values()

    for target in matched_targets:
        try:
            click.echo('Cancelling target {}.'.format(target.name))
            backend.cancel(target)
        except UnknownTargetError as exc:
            click.echo('Target {} could not be cancelled since it is unknown to the backend.'.format(exc))
        except UnsupportedOperationError:
            click.echo('Cancelling targets is not supported by this backend.')
            raise click.Abort()
