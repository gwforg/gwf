import click

from .. import Workflow
from ..backends import create_backend
from ..backends.exceptions import TargetError, UnsupportedOperationError
from ..core import CachedFilesystem, Graph, pass_context
from ..filtering import filter_names


def cancel_many(backend, targets):
    for target in targets:
        try:
            click.echo("Cancelling target {}".format(target.name), err=True)
            backend.cancel(target)
        except UnsupportedOperationError:
            click.echo("Cancelling targets is not supported by this backend", err=True)
            raise click.Abort()
        except TargetError:
            click.echo(
                f"Target {target.name} could not be cancelled "
                "(maybe not running or submitted?)"
            )


@click.command()
@click.argument("targets", nargs=-1)
@click.option(
    "-f", "--force", is_flag=True, default=False, help="Do not ask for confirmation."
)
@pass_context
def cancel(ctx, targets, force):
    """Cancel the specified targets."""

    if not force and not targets:
        click.confirm(
            "This will cancel all targets! Do you want to continue?", abort=True
        )

    workflow = Workflow.from_context(ctx)
    fs = CachedFilesystem()
    graph = Graph.from_targets(workflow.targets, fs)

    if targets:
        targets = filter_names(graph, targets)
    else:
        targets = list(graph)

    with create_backend(
        ctx.backend, working_dir=ctx.working_dir, config=ctx.config
    ) as backend:
        cancel_many(backend, targets)
