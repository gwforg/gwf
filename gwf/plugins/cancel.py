import click

from ..backends import backend_from_config
from ..backends.exceptions import UnsupportedOperationError, UnknownTargetError
from ..core import graph_from_config
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
@click.option('-f', '--force', is_flag=True, default=False, help='Do not ask for confirmation.')
@click.pass_obj
def cancel(obj, targets, force):
    """Cancel the specified targets."""
    graph = graph_from_config(obj)

    backend_cls = backend_from_config(obj)
    with backend_cls() as backend:
        if not targets and not force:
            if click.confirm('This will cancel all targets! Do you want to continue?', abort=True):
                cancel_many(backend, graph)
        else:
            matched_targets = filter_names(graph, targets)
            cancel_many(backend, matched_targets)
