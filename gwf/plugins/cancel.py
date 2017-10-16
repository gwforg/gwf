import click

from ..backends import backend_from_config
from ..core import graph_from_config
from ..exceptions import UnknownTargetError, UnsupportedOperationError
from ..filtering import filter_names


def cancel_many(backend, targets):
    for target in targets:
        try:
            click.echo('Cancelling target {}.'.format(target.name))
            backend.cancel(target)
        except UnknownTargetError as exc:
            click.echo('Target {} could not be cancelled since it is unknown to the backend.'.format(exc))
        except UnsupportedOperationError:
            click.echo('Cancelling targets is not supported by this backend.')
            raise click.Abort()


@click.command()
@click.argument('targets', nargs=-1)
@click.pass_obj
def cancel(obj, targets):
    """Cancel the specified targets."""
    graph = graph_from_config(obj)

    backend_cls = backend_from_config(obj)
    backend = backend_cls()

    matched_targets = filter_names(graph.targets.values(), targets)
    if not targets:
        if click.confirm('This will cancel all targets! Do you want to continue?', abort=True):
            cancel_many(backend, matched_targets)
    else:
        cancel_many(backend, matched_targets)
