"""
Implementation of the status command.
"""

import os
import colorama
import statusbar
from ..utils import dfs
from ..exceptions import TargetDoesNotExistError
from .base import Plugin


class StatusCommand(Plugin):

    name = "status"

    def __init__(self):
        self.ts = os.get_terminal_size()

    def _split_target_list(self, targets):
        should_run, submitted, running, completed = [], [], [], []
        for target in targets:
            if self.workflow.should_run(target):
                should_run.append(target)
            # FIXME: check for queue status here...how do I get the backend?
            else:
                completed.append(target)

        return should_run, submitted, running, completed

    def print_verbose(self, target_names):  # pragma: no cover
        columns = self.ts.columns
        status_string = " {{:.<{}}} {{:^10}}".format(columns - 13)

        for target_name in target_names:
            target = self.workflow.targets[target_name]
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = \
                self._split_target_list(dependencies)

            print(" {}".format(colorama.Style.BRIGHT + target_name + colorama.Style.NORMAL))
            print("_" * columns)
            for t in completed:
                print(status_string.format(
                    t.name,
                    colorama.Fore.GREEN + "DONE" + colorama.Fore.RESET
                ))
            for t in submitted:
                print(status_string.format(
                    t.name,
                    colorama.Fore.BLUE + "SUBMITTED" + colorama.Fore.RESET
                ))
            for t in running:
                print(status_string.format(
                    t.name,
                    colorama.Fore.YELLOW + "RUNNING" + colorama.Fore.RESET
                ))
            for t in should_run:
                print(status_string.format(
                    t.name,
                    colorama.Fore.RED + "SHOULD RUN" + colorama.Fore.RESET
                ))
            print("=" * columns)
            print()

    def print_progress(self, target_names):  # pragma: no cover
        table = statusbar.StatusTable()
        # FIXME: make status table here.
        # for target_name in target_names:
        #     target = self.workflow.targets[target_name]
        #     dependencies = dfs(target, self.workflow.dependencies)
        #     should_run, submitted, running, completed = self._split_target_list(dependencies)
        #     status_bar = make_status_bar(should_run, submitted, running, completed)
        #     print(status_string.format(Style.BRIGHT + target_name + Style.NORMAL, status_bar))
        print("\n".join(table.format_table()))

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            "status",
            "Command for getting the status of workflow targets.",
            self.on_run
        )

        subparser.add_argument("targets", metavar="TARGET", nargs="*",
                               help="Targets to show the status of (default all terminal targets)")
        subparser.add_argument("--verbose", action="store_true",
                               help="Output verbose status output")

    def on_run(self):
        targets = []
        if not self.config['targets']:
            targets = self.workflow.endpoints()
        else:
            for name in self.config['targets']:
                if name not in self.workflow.targets:
                    raise TargetDoesNotExistError(name)
                targets.append(self.workflow.targets[name])

        if self.config['verbose']:
            self.print_verbose(targets)
        else:
            self.print_progress(targets)
