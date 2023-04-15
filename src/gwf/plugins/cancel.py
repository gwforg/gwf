import click

from ..backends import Backend, Status
from ..backends.exceptions import UnsupportedOperationError
from ..core import CachedFilesystem, Graph
from ..filtering import filter_names


def cancel_many(backend, targets):
    for target in targets:
        try:
            click.echo("Cancelling target {}".format(target.name), err=True)
            if backend.status(target) != Status.UNKNOWN:
                backend.cancel(target)
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

    if not force and not targets:
        click.confirm(
            "This will cancel all targets! Do you want to continue?", abort=True
        )

    workflow = obj["workflow"]
    fs = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, fs)

    if targets:
        targets = filter_names(graph, targets)
    else:
        targets = list(graph)

    backend_cls = Backend.from_config(obj)
    with backend_cls() as backend:
        cancel_many(backend, targets)
