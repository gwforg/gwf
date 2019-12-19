import click

from ..backends import Backend
from ..backends.exceptions import TargetError, UnsupportedOperationError
from ..core import Graph
from ..filtering import filter_names


def cancel_many(backend, targets, ignore_unknown=False):
    for target in targets:
        try:
            click.echo("Cancelling target {}".format(target.name), err=True)
            backend.cancel(target)
        except TargetError:
            if ignore_unknown:
                continue

            click.echo(
                "Target {} could not be cancelled since it is unknown to the backend".format(
                    target.name
                ),
                err=True,
            )
        except UnsupportedOperationError:
            click.echo("Cancelling targets is not supported by this backend", err=True)
            raise click.Abort()


@click.command()
@click.argument("targets", nargs=-1)
@click.option(
    "-f", "--force", is_flag=True, default=False, help="Do not ask for confirmation."
)
@click.pass_obj
def cancel(obj, targets, force):
    """Cancel the specified targets."""
    graph = Graph.from_config(obj)
    backend_cls = Backend.from_config(obj)

    with backend_cls() as backend:
        if not targets:
            if not force:
                click.confirm(
                    "This will cancel all targets! Do you want to continue?", abort=True
                )
            cancel_many(backend, graph, ignore_unknown=True)
        else:
            matched_targets = filter_names(graph, targets)
            cancel_many(backend, matched_targets)
