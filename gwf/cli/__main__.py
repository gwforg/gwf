from __future__ import absolute_import, print_function
import sys
import colorama
from .arg_parsing import ArgumentDispatching
from .arg_parsing import SubCommand
from .status import StatusCommand

class RunCommand(SubCommand):
    def set_arguments(self, parser):
        parser.add_argument("targets", metavar = "TARGET", nargs = "*",
                            help = "Targets to run (default all terminal targets)")

    def handle(self, arguments):
        print("run: " + str(arguments))
        print(self.workflow)


def main():
    colorama.init()
    dispatcher = ArgumentDispatching()
    dispatcher.install_subcommand("run", "submit targets", RunCommand())
    dispatcher.install_subcommand("status", "check status of targets", StatusCommand())
    try:
        dispatcher.dispatch(sys.argv[1:])
    except Exception as e:
        print("Error: " + str(e))
