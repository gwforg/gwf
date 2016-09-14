from __future__ import absolute_import, print_function
import sys
from .arg_parsing import ArgumentDispatching
from .arg_parsing import SubCommand

class RunCommand(SubCommand):
    def set_arguments(self, parser):
        parser.add_argument("targets", metavar = "TARGET", nargs = "*",
                            help = "Targets to run (default all terminal targets)")

    def handle(self, arguments):
        print("run: " + str(arguments))
        print(self.workflow)

class StatusCommand(SubCommand):
    def set_arguments(self, parser):
        parser.add_argument("targets", metavar = "TARGET", nargs = "*",
                            help = "Targets to show status of (default all terminal targets)")

    def handle(self, arguments):
        print("status: " + str(arguments))
        print(self.workflow)

def main():
    dispatcher = ArgumentDispatching()
    dispatcher.install_subcommand("run", "submit targets", RunCommand())
    dispatcher.install_subcommand("status", "check status of targets", StatusCommand())
    dispatcher.dispatch(sys.argv[1:])
