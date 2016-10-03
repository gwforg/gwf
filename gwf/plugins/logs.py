import logging

from . import Plugin
from ..exceptions import TargetDoesNotExistError

logger = logging.getLogger(__name__)


class LogsCommand(Plugin):

    name = 'logs'
    help_text = (
        'Display logs for the latest run of a target. By default only '
        'standard output is shown. Supply the --stderr flag to show standard '
        'error instead.'
    )

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'logs',
            self.help_text,
            self.on_run
        )

        subparser.add_argument(
            '-e',
            '--stderr',
            action='store_true',
            help='Return stderr.',
        )

        subparser.add_argument(
            'target',
            help='Name of target.',
        )

    def on_run(self):
        target_name = self.config['target']
        if target_name not in self.workflow.targets:
            raise TargetDoesNotExistError(target_name)

        target = self.workflow.targets[target_name]
        log = self.backend.logs(target, stderr=self.config['stderr'])
        print(log.read(), end='')
