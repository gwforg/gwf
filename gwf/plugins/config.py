import logging
from collections import OrderedDict

from configargparse import DefaultConfigFileParser

from .. import LOCAL_CONFIG_FILE, USER_CONFIG_FILE
from ..exceptions import GWFError
from .base import Plugin

logger = logging.getLogger(__name__)


class ConfigCommand(Plugin):

    name = 'config'
    help_text = (
        'Set and unset configuration options. If a value is provided, the '
        'will be set to this value. If no value is given, the setting will be '
        'unset. To apply the change to the user configuration file, use the '
        '--user flag.'
    )

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'config',
            self.help_text,
            self.on_run
        )

        subparser.add_argument(
            '-u',
            '--user',
            action='store_true',
            help="Changes should apply to the user configuration."
        )

        subparser.add_argument(
            'option_name',
            metavar='option',
            help='Name of option.',
        )

        subparser.add_argument(
            'option_value',
            metavar='value',
            help='Value to set option to.',
            nargs='?',
            default='',
        )

    def _parse(self, parser, path, mode='r'):
        try:
            with open(path) as fileobj:
                return parser.parse(fileobj.readlines())
        except IOError as e:
            return OrderedDict()

    def on_run(self):
        option_name = self.config['option_name']
        option_value = self.config['option_value']
        user = self.config['user']

        config_path = USER_CONFIG_FILE if user else LOCAL_CONFIG_FILE

        parser = DefaultConfigFileParser()
        settings = self._parse(parser, config_path)

        if option_value:
            settings[option_name] = option_value
        else:
            if option_name not in settings:
                raise GWFError(
                    'Option could not be unset since it is not currently set.'
                )
            del settings[option_name]

        with open(config_path, 'w') as fileobj:
            fileobj.write(parser.serialize(settings))
