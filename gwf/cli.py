import logging
import platform
import sys

import colorama

from configargparse import ArgParser

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

    USER_CONFIG_FILE = '~/.gwfrc'
    LOCAL_CONFIG_FILE = '.gwfrc'

    def __init__(self):
        self.plugins_manager = ExtensionManager(group='gwf.plugins')
        self.backends_manager = ExtensionManager(group='gwf.backends')

        self.backend = None
        self.workflow = None

        self.parser = ArgParser(
            description=self.description,
            epilog=self.epilog,
            default_config_files=[
                self.USER_CONFIG_FILE, self.LOCAL_CONFIG_FILE
            ],
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
            default='slurm',
        )

        self.parser.add_argument(
            '-c',
            '--config',
            is_config_file=True,
            help='path to a specific configuration file.'
        )

        self.parser.add_argument(
            '-v',
            '--verbosity',
            help='set verbosity level.',
            default='warning',
        )

        # Prepare for sub-commands
        self.subparsers = \
            self.parser.add_subparsers(title="commands")

    def configure_logging(self):
        logging.basicConfig(
            level={
                'warning': logging.WARNING,
                'info': logging.INFO,
                'debug': logging.DEBUG,
            }[self.args.verbosity],
            format='%(levelname)-6s|  %(message)s',
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
        logger.debug('Config:\n\n%s', self.parser.format_values())

        self.load_workflow()

        self.backend = self.backends_manager.exts[self.args.backend]
        self.backend.configure(
            workflow=self.prepared_workflow,
            config=self.args,
        )

        self.plugins_manager.configure_extensions(
            workflow=self.prepared_workflow,
            backend=self.backend,
            config=self.args,
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
