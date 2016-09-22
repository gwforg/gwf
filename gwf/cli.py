import argparse
import logging
import sys

import colorama
from pkg_resources import iter_entry_points

from .core import PreparedWorkflow
from .exceptions import GWFError
from .utils import import_object

logger = logging.getLogger(__name__)


class App:

    description = "A flexible, pragmatic workflow tool."

    def __init__(self):
        self.config = {'defaults': {'cores': 16}, 'backend': 'testing'}

        self.plugin_classes = self.load_extensions('gwf.plugins')
        self.backend_classes = self.load_extensions('gwf.backends')

        self.plugins = {}
        self.backend = None

        self.parser = argparse.ArgumentParser(
            description=self.description,
        )

        # Set global options here. Options for sub-commands will be set by the
        # sub-commands we dispatch to.
        self.parser.add_argument(
            "-f",
            "--file",
            help="workflow file/object (default: workflow.py).",
            default="workflow.py:gwf"
        )

        self.parser.add_argument(
            '-b',
            '--backend',
            help='backend used to run the workflow.',
            choices=self.backend_classes.keys(),
            default=self.config['backend'],
        )

        self.parser.add_argument(
            '-v',
            '--verbosity',
            help='set verbosity level (default: 1).',
            default=1,
            action='count',
        )

        # Prepare for sub-commands
        self.subparsers = \
            self.parser.add_subparsers(title="Commands")

    def load_extensions(self, group):
        exts = {}
        for entry_point in iter_entry_points(group=group, name=None):
            ext_cls = entry_point.load()
            logger.debug('Loaded extension: %s.', ext_cls.name)
            exts[ext_cls.name] = ext_cls
        return exts

    def init_plugins(self):
        for plugin_name, plugin_cls in self.plugin_classes.items():
            logger.debug('Initializing plugin: %s.', plugin_name)

            plugin = plugin_cls()
            plugin.setup_argument_parser(self.parser, self.subparsers)
            self.plugins[plugin.name] = plugin

    def configure_plugins(self):
        for plugin_name, plugin in self.plugins.items():
            logger.debug('Configuring plugin: %s.', plugin_name)

            plugin.configure(
                workflow=self.prepared_workflow,
                backend=self.backend,
                config=self.config,
                args=self.args,
            )

    def run(self, argv):
        """Parse command line arguments in argv and dispatch to subcommand."""

        # Initialize plugins.
        self.init_plugins()

        # Parse arguments.
        self.args = self.parser.parse_args(argv)

        # Configure logging.
        logging.basicConfig(
            level={
                1: logging.WARNING,
                2: logging.INFO,
                3: logging.DEBUG,
            }[self.args.verbosity],
            format='[%(levelname)s] %(message)s',
        )

        # Read and prepare the workflow.
        logging.debug('Loading workflow from: %s.', self.args.file)
        workflow = import_object(self.args.file)
        self.prepared_workflow = PreparedWorkflow(workflow=workflow)

        # Initialize and configure backend.
        backend_cls = self.backend_classes[self.args.backend]

        logger.debug('Initializing backend: %s.', backend_cls.name)
        self.backend = backend_cls()

        logger.debug('Configuring backend: %s.', backend_cls.name)
        self.backend.configure(
            workflow=self.prepared_workflow,
            config=self.config,
            args=self.args,
        )

        # Initialize and configure all plugins.
        self.configure_plugins()

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
        print("[Error] {}".format(str(e)))
        sys.exit(1)
