import json
from collections import ChainMap


CONFIG_DEFAULTS = {
    'verbose': 'info',
    'backend': 'local',
    'no_color': False,
}


class FileConfig(dict):

    def __init__(self, path, initial):
        super().__init__(initial)
        self.path = path

    def dump(self):
        """Dump the configuration to disk."""
        with open(self.path, 'w') as config_file:
            json.dump(self, config_file, indent=4, sort_keys=True)

    @classmethod
    def load(cls, path):
        """Load configuration from a file.

        Reads configuration from `file` and returns a :class:`Config` instance
        with the configuration. The `defaults` will be merged into the
        configuration.

        :param path str: Path to the configuration file.
        :param defaults dict: A set of defaults to merge into the configuration.
        """
        try:
            with open(path) as config_file:
                file_config = json.load(config_file)
        except FileNotFoundError:
            file_config = {}
        return cls(path=path, initial=file_config)


file_config = FileConfig.load('.gwfconf.json')
config = ChainMap(file_config, CONFIG_DEFAULTS)
