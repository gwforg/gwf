from __future__ import absolute_import, print_function
from .arg_parsing import SubCommand
from ..utils import dfs
from ..exceptions import TargetDoesNotExistsError
from colorama import Fore, Back, Style
import os
from math import ceil

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
        columns = self.ts.columns
        status_string = " {{:.<{}}} {{:^10}}".format(columns - 13)

        for target_name in target_names:
            target = self.workflow.targets[target_name]
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = self._split_target_list(dependencies)

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


    def print_progress(self, target_names):
        columns = self.ts.columns
        name_width = status_width = int((columns) / 2)
        status_string = " {{:.<{}}} {{:^{}}}".format(name_width, status_width)

        def get_width(left, n, k):
            if k == 0: return 0
            return int(ceil(left * n / k))

        def make_status_bar(should_run, submitted, running, completed):
            n_should_run = len(should_run)
            n_submitted = len(submitted)
            n_running = len(running)
            n_completed = len(completed)
            n_total = n_should_run + n_submitted + n_running + n_completed

            # I am using two characters for brackets, so I have status_width - 2
            # characters to fill out. The less complete a job is, the more important
            # it is to show it.
            n = left = int(status_width - 2)
            should_run_width = get_width(left, n_should_run, n_should_run + n_submitted + n_running + n_completed)
            left -= should_run_width
            submitted_width = get_width(left, n_submitted, n_submitted + n_running + n_completed)
            left -= submitted_width
            running_width = get_width(left, n_running, n_running + n_completed)
            left -= running_width
            completed_width = left

            completed_bar = Fore.GREEN + ("#" * completed_width) + Fore.RESET
            running_bar = Fore.YELLOW + ("#" * running_width) + Fore.RESET
            submitted_bar = Fore.BLUE + ("." * submitted_width) + Fore.RESET
            should_run_bar = Fore.RED + ("." * should_run_width) + Fore.RESET

            return "[{}]".format(completed_bar + running_bar + submitted_bar + should_run_bar)


        for target_name in target_names:
            target = self.workflow.targets[target_name]
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = self._split_target_list(dependencies)
            status_bar = make_status_bar(should_run, submitted, running, completed)
            print(status_string.format(Style.BRIGHT + target_name + Style.NORMAL, status_bar))


    def handle(self, arguments):
        target_names = arguments.targets
        if len(target_names) == 0:
            target_names = [target.name for target in self.workflow.endpoints()]

        # Check targets are in workflow
        for target_name in target_names:
            if target_name not in self.workflow.targets:
                raise TargetDoesNotExistsError(target_name)

        if arguments.verbose:
            self.print_verbose(target_names)
        else:
            self.print_progress(target_names)
