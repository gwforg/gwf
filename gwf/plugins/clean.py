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

        subparser.add_argument(
            '-e',
            '--not-endpoints',
            action='store_true',
            help='Only delete output files from targets that are not endpoints.'
        )

    def on_clean(self):
        workflow = self.get_prepared_workflow()
        backend = self.get_active_backend()

        targets = set()
        if not self.config['targets']:
            targets.update(workflow.targets.values())
        else:
            for name in self.config['targets']:
                if name not in workflow.targets:
                    raise TargetDoesNotExistError(name)
                targets.add(workflow.targets[name])

        if self.config['not_endpoints']:
            targets -= workflow.endpoints()

        for target in targets:
            if not self.config['only_failed'] or backend.failed(target):
                logger.info('Cleaning up: %s.', target.name)
                for path in target.outputs:
                    logging.debug(
                        'Deleting output file "%s" from target "%s".',
                        path,
                        target.name
                    )

                    _delete_file(path)
