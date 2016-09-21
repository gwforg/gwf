from __future__ import absolute_import, print_function
from ..utils import cache, import_object
from ..core import PreparedWorkflow
import argparse


class SubCommand(object):
    """Super-class for sub-commands."""

    # Class variable used to get the active workflow.
    _workflow = None
    @property
    @cache
    def workflow(cls):
        if cls._workflow is None: # pragma: no cover
            sys.exit("Workflow not specified") # FIXME: exception, but it really should never happen since we have a default
        return PreparedWorkflow(import_object(cls._workflow))

    def set_arguments(self, subparser):
        "Provide the parser options for the sub-command."

    def handle(self, arguments):
        "Callback when this sub-command is invoked."

class ArgumentDispatching(object):
    """Class responsible for setting up the command line arguments and dispatching
    to sub-commands when invoked."""

    def __init__(self):
        self.parser = argparse.ArgumentParser(
            description="Grid WorkFlow",
        )

        # Set global options here. Options for sub-commands will be set by the
        # sub-commands we dispatch to.
        self.parser.add_argument("-w", "--workflow",
                                 help = "workflow file/object (default 'workflow.py')",
                                 default = "workflow.py")

        # Prepare for sub-commands
        self.subparsers = self.parser.add_subparsers(title = "commands", description = "Subcommands to execute")

    def install_subcommand(self, name, description, handler):
        """Install a sub-command with "name" and with functionality implemented
        in "handler"."""
        subparser = self.subparsers.add_parser(name, help = description)
        handler.set_arguments(subparser)
        subparser.set_defaults(func = handler.handle)

    def dispatch(self, argv):
        """Parse command line arguments in argv and dispatch to subcommand."""
        args = self.parser.parse_args(argv)
        SubCommand._workflow = args.workflow
        if hasattr(args, "func"):
            args.func(args)
