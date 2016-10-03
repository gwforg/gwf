import logging

from . import Plugin
from ..exceptions import TargetDoesNotExistError

logger = logging.getLogger(__name__)


class CleanCommand(Plugin):

    name = 'clean'

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'clean',
            'Clean up output files from one or more targets.',
            self.on_clean
        )

        subparser.add_argument(
            "targets",
            metavar="TARGET",
            nargs="*",
            help="Targets to clean up (default: all targets)"
        )

    def on_clean(self):
        targets = []
        if not self.config['targets']:
            targets = self.workflow.targets.values()
        else:
            for name in self.config['targets']:
                if name not in self.workflow.targets:
                    raise TargetDoesNotExistError(name)
                targets.append(self.workflow.targets[name])

        for target in targets:
            logger.info('Cleaning up: %s.', target.name)
            target.clean()
