import os
import textwrap

import statusbar

from . import Plugin
from ..exceptions import TargetDoesNotExistError
from ..utils import dfs


class StatusCommand(Plugin):
    """
    Show the status of targets.

    By default, shows a progress bar for each endpoint in the workflow.
    If one or more target names are supplied, progress bars are shown
    for these targets.

    A progress bar represents the target and its dependencies, and
    shows how many of the dependencies either should run (red, .),
    are submitted (yellow, -), are running (blue, +) or are
    completed (green, #).
    """

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
                elif self.backend.completed(target) or self.backend.failed(target):
                    completed.append(target)
            else:
                completed.append(target)

        return should_run, submitted, running, completed

    def print_progress(self, targets):  # pragma: no cover
        table = statusbar.StatusTable()
        for target in targets:
            dependencies = dfs(target, self.workflow.dependencies)
            should_run, submitted, running, completed = self._split_target_list(
                dependencies)
            status_bar = table.add_status_line(target.name)
            status_bar.add_progress(len(completed), '#', color='green')
            status_bar.add_progress(len(running), '+', color='blue')
            status_bar.add_progress(len(submitted), '-', color='yellow')
            status_bar.add_progress(len(should_run), '.', color='red')
        print('\n'.join(table.format_table()))

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'status',
            textwrap.dedent(StatusCommand.__doc__),
            self.on_run)

        subparser.add_argument('targets', metavar='TARGET', nargs='*',
                               help='Targets to show the status of (default: all terminal targets)')

    def on_run(self):
        self.workflow = self.get_prepared_workflow()
        self.backend = self.get_active_backend()

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
        self.print_progress(targets)
