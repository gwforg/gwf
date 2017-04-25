import logging
from argparse import RawDescriptionHelpFormatter

logger = logging.getLogger(__name__)


class Extension:

    @staticmethod
    def setup_subparser(subparsers, name, description, handler):
        """Helper method for setting up subparsers.

        :param subparsers:
            A subparsers object as returned by
            :func:`ArgumentParser.add_subparsers`.
        :param name:
            Name of the subcommand.
        :param description:
            A short description of the subcommand.
        :param handler:
            The callable that will be called when the subcommand is used.

        This method is most often used in :func:`setup_argument_parser` to set
        up subcommands. For example, to register a subcommand called ``foo``
        that calls the method ``on_foo``::

            class MyPlugin(Plugin):
                name = 'myplugin'

                def setup_argument_parser(self, parser, subparsers):
                    self.setup_subparser(
                        subparsers,
                        'foo',
                        'a command that prints a greeting.',
                        self.on_foo,
                    )

                def on_foo(self):
                    print('Hello!')

        When the plugin is installed and registered properly, the user will now
        be able to run ``gwf foo`` which will print ``Hello!`` and exit.
        """
        subparser = subparsers.add_parser(
            name,
            description=description,
            formatter_class=RawDescriptionHelpFormatter)
        subparser.set_defaults(func=handler)
        return subparser
