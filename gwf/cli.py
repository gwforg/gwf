import argparse
import logging
import platform
import sys

import colorama

from .conf import settings
from .core import PreparedWorkflow
from .exceptions import GWFError
from .ext import ExtensionManager
from .utils import get_gwf_version, import_object

logger = logging.getLogger(__name__)


class App:

    description = "A flexible, pragmatic workflow tool."

    epilog = (
        'gwf is a flexible, pragmatic workflow tool for building and running '
        'large, scientific workflows developed at the Bioinformatics Research '
        'Centre (BiRC), Aarhus University.'
    )

    def __init__(self):
        self.plugins_manager = ExtensionManager(group='gwf.plugins')
        self.backends_manager = ExtensionManager(group='gwf.backends')

        self.backend = None
        self.workflow = None

        self.parser = argparse.ArgumentParser(
            description=self.description,
            epilog=self.epilog,
        )

        # Set global options here. Options for sub-commands will be set by the
        # sub-commands we dispatch to.
        self.parser.add_argument(
            "-f",
            "--file",
            help="workflow file/object (default: workflow.py:gwf).",
            default="workflow.py:gwf"
        )

        self.parser.add_argument(
            '-b',
            '--backend',
            help='backend used to run the workflow.',
            choices=self.backends_manager.exts.keys(),
            default=settings['backend'],
        )

        self.parser.add_argument(
            '-v',
            '--verbosity',
            help='increase verbosity level.',
            default=1,
            action='count',
        )

        # Prepare for sub-commands
        self.subparsers = \
            self.parser.add_subparsers(title="commands")

    def configure_logging(self):
        logging.basicConfig(
            level={
                1: logging.WARNING,
                2: logging.INFO,
                3: logging.DEBUG,
            }[self.args.verbosity],
            format='[%(levelname)s] %(message)s',
        )

    def load_workflow(self):
        logger.debug('Loading workflow from: %s.', self.args.file)
        workflow = import_object(self.args.file)
        self.prepared_workflow = PreparedWorkflow(workflow=workflow)

    def run(self, argv):
        self.plugins_manager.load_extensions()
        self.backends_manager.load_extensions()

        self.plugins_manager.setup_argument_parsers(
            self.parser, self.subparsers)
        self.backends_manager.setup_argument_parsers(
            self.parser, self.subparsers)

        self.args = self.parser.parse_args(argv)

        self.configure_logging()

        logger.debug('Platform: %s.', platform.platform())
        logger.debug('GWF version: %s.', get_gwf_version())
        logger.debug('Python version: %s.', platform.python_version())
        logger.debug('Node: %s.', platform.node())

        self.load_workflow()

        self.backend = self.backends_manager.exts[self.args.backend]
        self.backend.configure(
            workflow=self.prepared_workflow,
            config=settings,
            args=self.args,
        )

        self.plugins_manager.configure_extensions(
            workflow=self.prepared_workflow,
            backend=self.backend,
            config=settings,
            args=self.args,
        )

        # Dispatch to subcommand.
        if hasattr(self.args, "func"):
            self.args.func()
        else:
            self.parser.print_help()
            sys.exit(1)


def main():
    colorama.init()

    app = App()
    try:
        app.run(sys.argv[1:])
    except GWFError as e:
        print("[ERROR] {}".format(str(e)))
        sys.exit(1)
