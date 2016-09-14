from __future__ import absolute_import, print_function
from .arg_parsing import SubCommand
from ..utils import dfs
from colorama import Fore, Back, Style
import os

class StatusCommand(SubCommand):
    def _split_target_list(self, targets):
        should_run, submitted, running, completed = [], [], [], []
        for target in targets:
            if self.workflow.should_run(target):
                should_run.append(target)
            #FIXME: check for queue status here...how do I get the backend?
            else:
                completed.append(target)

        return should_run, submitted, running, completed

    def __init__(self):
        self.ts = ts = os.get_terminal_size()

    def set_arguments(self, parser):
        parser.add_argument("targets", metavar = "TARGET", nargs = "*",
                            help = "Targets to show the status of (default all terminal targets)")
        parser.add_argument("--verbose", action = "store_true",
                            help = "Output verbose status output")

    def print_verbose(self, target_names):
        for target_name in target_names:
            target = self.workflow.targets[target_name]
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = self._split_target_list(dependencies)

            columns = self.ts.columns
            status_string = " {{:.<{}}} {{:^10}}".format(columns - 13)

            #print()
            #print("_" * columns)
            print(" {}".format(Style.BRIGHT + target_name + Style.NORMAL))
            print("_" * columns)
            for t in completed:
                print(status_string.format(t.name, Fore.GREEN + "DONE" + Fore.RESET))
            for t in submitted:
                print(status_string.format(t.name, Fore.BLUE + "SUBMITTED" + Fore.RESET))
            for t in running:
                print(status_string.format(t.name, Fore.YELLOW + "RUNNING" + Fore.RESET))
            for t in should_run:
                print(status_string.format(t.name, Fore.RED + "SHOULD RUN" + Fore.RESET))
            print("=" * columns)
            print()


    def handle(self, arguments):
        target_names = arguments.targets
        if len(target_names) == 0:
            target_names = [target.name for target in self.workflow.endpoints()]

        # FIXME: check targets are in workflow

        if arguments.verbose:
            self.print_verbose(target_names)
