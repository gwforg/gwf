import sys
import colorama
from .parsing import Dispatcher
from .parsing import Command
from .status import StatusCommand


class RunCommand(Command):
    def set_arguments(self, parser):
        parser.add_argument("targets", metavar="TARGET", nargs="*",
                            help="Targets to run (default all terminal targets)")

    def handle(self, arguments):
        print("run: " + str(arguments))
        print(self.workflow)


def main():
    colorama.init()
    dispatcher = Dispatcher()
    dispatcher.install_subcommand("run", "submit targets", RunCommand())
    dispatcher.install_subcommand("status", "check status of targets", StatusCommand())
    try:
        dispatcher.dispatch(sys.argv[1:])
    except Exception as e:
        print("Error: " + str(e))
