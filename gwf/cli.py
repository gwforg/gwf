import logging
import os
import os.path
import platform
import sys

import colorama

from configargparse import ArgParser

from .core import PreparedWorkflow
from .exceptions import GWFError
from .ext import ExtensionManager
from .utils import get_gwf_version, import_object, merge

logger = logging.getLogger(__name__)


USER_CONFIG_FILE = '~/.gwfrc'
LOCAL_CONFIG_FILE = '.gwfrc'


class App:

    description = "A flexible, pragmatic workflow tool."

    VERBOSITY_LEVELS = {
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG,
    }

    def __init__(self, plugins_manager, backends_manager, config_files):
        self.plugins_manager = plugins_manager
        self.backends_manager = backends_manager

        self.active_backend = None
        self.workflow = None

        self.parser = ArgParser(
            description=self.description,
            add_config_file_help=False,
            add_env_var_help=False,
            auto_env_var_prefix='GWF_',
            default_config_files=config_files,
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
            choices=self.VERBOSITY_LEVELS.keys()
        )

        # Prepare for subcommands
        self.subparsers = \
            self.parser.add_subparsers(title="commands")

        # Load plugins and register their arguments and subcommands.
        self.plugins_manager.load_extensions()
        self.backends_manager.load_extensions()

        self.plugins_manager.setup_argument_parsers(
            self.parser, self.subparsers
        )
        self.backends_manager.setup_argument_parsers(
            self.parser, self.subparsers
        )

    def _configure_logging(self):
        logging.basicConfig(
            level=self.VERBOSITY_LEVELS[self.config['verbosity']],
            format='%(levelname)-6s|  %(message)s',
        )

    def run(self, argv):
        self.config = vars(self.parser.parse_args(argv))

        # Add path of workflow file to python path to make it possible to load
        # modules placed in the directory directly.
        workflow_dir = os.path.dirname(os.path.abspath(self.config['file']))
        sys.path.insert(0, workflow_dir)

        # Make sure that a .gwf directory exists in the workflow directory.
        try:
            os.mkdir(os.path.join(workflow_dir, '.gwf'))
        except IOError:
            raise GWFError(
                'Could not create directory in workflow directory "{}".'.format(
                    workflow_dir
                )
            )

        # If a subcommand is being called, the handler will be the function to
        # call when all loading is done.
        handler = self.config.pop('func', None)

        self._configure_logging()

        logger.debug('Platform: %s.', platform.platform())
        logger.debug('GWF version: %s.', get_gwf_version())
        logger.debug('Python version: %s.', platform.python_version())
        logger.debug('Node: %s.', platform.node())
        logger.debug('Python path: %s.', sys.path)

        backend_name = self.config['backend']
        logger.debug('Setting active backend: %s.', backend_name)
        self.active_backend = self.backends_manager.exts[backend_name]

        logger.debug('Loading workflow from: %s.', self.config['file'])
        workflow = import_object(self.config['file'])

        # Update configuration with defaults from backend and workflow.
        self.config = merge(
            self.active_backend.option_defaults,
            self.config,
            workflow.defaults,
        )

        self.prepared_workflow = PreparedWorkflow(
            workflow=workflow,
            config=self.config,
            backend=self.active_backend,
        )

        self.active_backend.configure(
            workflow=self.prepared_workflow,
            config=self.config,
        )

        self.plugins_manager.configure_extensions(
            workflow=self.prepared_workflow,
            backend=self.active_backend,
            config=self.config,
        )

        # Dispatch to subcommand.
        if handler is not None:
            handler()
        else:
            self.parser.print_help()

    def exit(self):
        logger.debug('Exiting...')
        self.plugins_manager.close_extensions()
        self.active_backend.close()


def main():
    colorama.init()

    config_files = [USER_CONFIG_FILE, LOCAL_CONFIG_FILE]

    plugins_manager = ExtensionManager(group='gwf.plugins')
    backends_manager = ExtensionManager(group='gwf.backends')

    app = App(
        plugins_manager=plugins_manager,
        backends_manager=backends_manager,
        config_files=config_files,
    )

    try:
        app.run(sys.argv[1:])
    except GWFError as e:
        logger.error(str(e))
    finally:
        app.exit()
