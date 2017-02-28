import os.path
import logging
import platform
import sys

from configargparse import ArgParser
from argparse import RawDescriptionHelpFormatter

from . import LOCAL_CONFIG_FILE, USER_CONFIG_FILE, __version__
from .core import PreparedWorkflow
from .exceptions import GWFError
from .ext import ExtensionManager
from .utils import cache, import_object, merge

logger = logging.getLogger(__name__)


class App:

    description = "A flexible, pragmatic workflow tool."

    VERBOSITY_LEVELS = {
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG,
    }

    BASIC_FORMAT = '%(levelname)-8s|  %(message)s'
    ADVANCED_FORMAT = '%(levelname)s:%(name)s:%(lineno)d| %(message)s'

    LOGGING_FORMATS = {
        'warning': BASIC_FORMAT,
        'info': BASIC_FORMAT,
        'debug': ADVANCED_FORMAT
    }

    def __init__(self, plugins_manager, backends_manager, config_files):
        self.plugins_manager = plugins_manager
        self.backends_manager = backends_manager

        self.active_backend = None
        self.workflow = None
        self.config = {}

        self.parser = ArgParser(
            description=self.description,
            add_config_file_help=False,
            add_env_var_help=False,
            auto_env_var_prefix='GWF_',
            default_config_files=config_files,
            formatter_class=RawDescriptionHelpFormatter,
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

        self.parser.add_argument(
            '--version',
            action='version',
            version=__version__
        )

        # Prepare for subcommands
        self.subparsers = self.parser.add_subparsers(title="commands")

        # Load plugins and register their arguments and subcommands.
        self.plugins_manager.load_extensions()
        self.backends_manager.load_extensions()

        self.plugins_manager.setup_argument_parsers(
            self.parser, self.subparsers
        )
        self.backends_manager.setup_argument_parsers(
            self.parser, self.subparsers
        )

    def _configure_logging(self, verbosity):
        logging.basicConfig(
            level=self.VERBOSITY_LEVELS[verbosity],
            format=self.LOGGING_FORMATS[verbosity],
        )

    def run(self, argv):
        self.config = vars(self.parser.parse_args(argv))
        self._configure_logging(self.config['verbosity'])

        # If a subcommand is being called, the handler will be the function to
        # call when all loading is done.
        handler = self.config.pop('func', None)

        logger.debug('Platform: %s.', platform.platform())
        logger.debug('GWF version: %s.', __version__)
        logger.debug('Python version: %s.', platform.python_version())
        logger.debug('Node: %s.', platform.node())
        logger.debug('Python path: %s.', sys.path)
        logger.debug(
            'Loaded plugins: %s.',
            ', '.join(self.plugins_manager.exts.keys())
        )
        logger.debug(
            'Loaded backends: %s.',
            ', '.join(self.backends_manager.exts.keys())
        )

        backend_name = self.config['backend']
        logger.debug('Setting active backend: %s.', backend_name)
        self.active_backend = self.backends_manager.exts[backend_name]

        logger.debug('Loading workflow from: %s.', self.config['file'])
        workflow = import_object(self.config['file'])

        # Add path of workflow file to python path to make it possible to load
        # modules placed in the directory directly.
        sys.path.insert(1, workflow.working_dir)

        # Ensure that a .gwf directory exists in the workflow directory.
        ensure_dir(os.path.join(workflow.working_dir, '.gwf'))

        # Update configuration with defaults from backend and workflow.
        self.config = merge(
            self.active_backend.option_defaults,
            self.config,
            workflow.defaults,
        )
        logger.debug('Merged configuration: %r.', self.config)

        @cache
        def get_prepared_workflow():
            return PreparedWorkflow(
                targets=workflow.targets,
                working_dir=workflow_dir,
                supported_options=self.active_backend.supported_options,
                config=self.config,
            )

        @cache
        def get_active_backend():
            self.active_backend.configure(
                working_dir=workflow_dir,
                config=self.config,
            )
            return self.active_backend

        self.plugins_manager.configure_extensions(
            get_prepared_workflow=get_prepared_workflow,
            get_active_backend=get_active_backend,
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
        if self.active_backend:
            self.active_backend.close()


def main():
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
