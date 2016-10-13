"""
Implementation of the status command.
"""

import os

import colorama

import statusbar

from . import Plugin
from ..exceptions import TargetDoesNotExistError
from ..utils import dfs


class StatusCommand(Plugin):

    def configure(self, *args, **kwargs):
        super().configure(*args, **kwargs)
        self.ts = os.get_terminal_size()

    def _split_target_list(self, targets):
        should_run, submitted, running, completed = [], [], [], []
        for target in targets:
            if self.workflow.should_run(target):
                should_run.append(target)
                if self.backend.submitted(target):
                    submitted.append(target)
                elif self.backend.running(target):
                    running.append(target)
            else:
                completed.append(target)

        return should_run, submitted, running, completed

    def print_verbose(self, targets):  # pragma: no cover
        columns = self.ts.columns
        status_string = " {{:.<{}}} {{:^10}}".format(columns - 13)

        for target in targets:
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = \
                self._split_target_list(dependencies)

            print(" {}".format(colorama.Style.BRIGHT +
                               target.name + colorama.Style.NORMAL))
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

    def print_progress(self, targets):  # pragma: no cover
        table = statusbar.StatusTable()
        for target in targets:
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = self._split_target_list(
                dependencies)
            status_bar = table.add_status_line(target.name + ": ")
            status_bar.add_progress(len(completed), "#", color="green")
            status_bar.add_progress(len(running), "#", color="blue")
            status_bar.add_progress(len(submitted), "-", color="yellow")
            status_bar.add_progress(len(should_run), ".", color="red")
        print("\n".join(table.format_table()))

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            "status",
            "Command for getting the status of workflow targets.",
            self.on_run
        )

        subparser.add_argument("targets", metavar="TARGET", nargs="*",
                               help="Targets to show the status of (default: all terminal targets)")
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

        # disabled coverage test on these since we don't actually test
        # the printing code...
        if self.config['verbose']:  # pragma: no cover
            self.print_verbose(targets)
        else:                      # pragma: no cover
            self.print_progress(targets)
