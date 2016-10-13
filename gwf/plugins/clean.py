import logging
import os

from . import Plugin
from ..exceptions import TargetDoesNotExistError

logger = logging.getLogger(__name__)


def _delete_file(path):
    try:
        os.remove(path)
    except OSError:
        pass


class CleanCommand(Plugin):

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

        subparser.add_argument(
            '-f',
            '--only-failed',
            action='store_true',
            help='Only delete output files from failed targets.',
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
            if not self.config['only_failed'] or self.backend.failed(target):
                logger.info('Cleaning up: %s.', target.name)
                for path in target.outputs:
                    logging.debug(
                        'Deleting output file "%s" from target "%s".',
                        path,
                        target.name
                    )

                    _delete_file(path)
