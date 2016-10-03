from . import Plugin
from ..exceptions import TargetDoesNotExistError


class RunCommand(Plugin):

    name = 'run'

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'run',
            'Command for running a workflow.',
            self.on_run
        )

        subparser.add_argument(
            "targets",
            metavar="TARGET",
            nargs="*",
            help="Targets to run (default: all terminal targets)"
        )

    def on_run(self):
        targets = []
        if not self.config['targets']:
            targets = self.workflow.endpoints()
        else:
            for name in self.config['targets']:
                if name not in self.workflow.targets:
                    raise TargetDoesNotExistError(name)
                targets.append(self.workflow.targets[name])
        self.backend.schedule_many(targets)
