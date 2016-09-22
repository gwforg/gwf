from ..ext import Plugin


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

    def configure(self, workflow, backend, config, args):
        self.workflow = workflow
        self.backend = backend
        self.config = config
        self.args = args

    def on_run(self):
        if not self.args.targets:
            targets = self.workflow.endpoints()
        else:
            targets = [self.workflow.targets[target]
                       for target in self.args.targets]

        self.backend.schedule_many(targets)
