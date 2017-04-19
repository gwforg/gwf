from . import Plugin
from ..core import schedule_many
from ..exceptions import TargetDoesNotExistError


class RunCommand(Plugin):

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
        workflow = self.get_graph()
        backend = self.get_active_backend()

        targets = []
        if not self.config['targets']:
            targets = workflow.endpoints()
        else:
            for name in self.config['targets']:
                if name not in workflow.targets:
                    raise TargetDoesNotExistError(name)
                targets.append(workflow.targets[name])
        schedule_many(workflow, backend, targets)
