from __future__ import absolute_import, print_function
from ..utils import cache, import_object
from ..core import PreparedWorkflow
import argparse
import sys

## Set default sub-command. Taken from https://bitbucket.org/ruamel/std.argparse
def set_default_subparser(self, name, args=None):
    """default subparser selection. Call after setup, just before parse_args()
    name: is the name of the subparser to call by default
    args: if set is the argument list handed to parse_args()

    , tested with 2.7, 3.2, 3.3, 3.4
    it works with 2.6 assuming argparse is installed
    """
    subparser_found = False
    for arg in sys.argv[1:]:
        if arg in ['-h', '--help']:  # global help if no subparser
            break
    else:
        for x in self._subparsers._actions:
            if not isinstance(x, argparse._SubParsersAction):
                continue
            for sp_name in x._name_parser_map.keys():
                if sp_name in sys.argv[1:]:
                    subparser_found = True
        if not subparser_found:
            # insert default in first position, this implies no
            # global options without a sub_parsers specified
            if args is None:
                sys.argv.insert(1, name)
            else:
                args.insert(0, name)

argparse.ArgumentParser.set_default_subparser = set_default_subparser


class SubCommand(object):

    # Class variable used to get the active workflow.
    _workflow = None
    @property
    @cache
    def workflow(cls):
        if cls._workflow is None:
            sys.exit("Workflow not specified") # FIXME: exception
        return PreparedWorkflow(import_object(cls._workflow))

    """Super-class for sub-commands."""
    def __init__(self):
        pass

    def set_arguments(self, subparser):
        pass

    def handle(self):
        pass

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
        self.subparsers = self.parser.add_subparsers(title = "commands",
                                                     description = "Subcommands to execute. The default is 'run'.")

    def install_subcommand(self, name, description, handler):
        """Install a sub-command with "name" and with functionality implemented
        in "handler"."""
        subparser = self.subparsers.add_parser(name, help = description)
        handler.set_arguments(subparser)
        subparser.set_defaults(func = handler.handle)

    def set_default_command(self, name):
        self.parser.set_default_subparser(name)

    def dispatch(self, argv):
        """Parse command line arguments in argv and dispatch to subcommand."""
        args = self.parser.parse_args(argv)
        SubCommand._workflow = args.workflow
        if hasattr(args, "func"):
            args.func(args)
