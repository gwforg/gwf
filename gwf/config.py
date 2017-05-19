import json
from collections import ChainMap

from .exceptions import GWFError


CONFIG_DEFAULTS = {
    'verbose': 'info',
    'backend': 'local',
    'no_color': False,
}


class Config:

    def __init__(self, path, config):
        self.path = path
        self._config = config

    @classmethod
    def load(cls, path, defaults):
        try:
            with open(path) as config_file:
                file_config = json.load(config_file)
        except FileNotFoundError:
            file_config = {}

        if defaults is None:
            defaults = {}
        return cls(path=path, config=ChainMap(file_config, defaults))

    def get_config(self):
        return dict(self._config)

    def get_file_config(self):
        return Config(path=self.path, config=self._config.maps[0])

    def dump(self):
        with open(self.path, 'w') as config_file:
            json.dump(dict(self._config), config_file, indent=4, sort_keys=True)

    def set(self, key, new_value):
        parts = key.split('.')
        config = self._config
        for subkey in parts[:-1]:
            if subkey not in config:
                config[subkey] = {}
            config = config[subkey]
        config[parts[-1]] = new_value

    def unset(self, key):
        parts = key.split('.')
        config = self._config
        try:
            for subkey in parts[:-1]:
                config = config[subkey]
            del config[parts[-1]]
        except KeyError:
            raise GWFError('Key does not exist.')

    def get(self, key, default=None):
        parts = key.split('.')
        try:
            value = self._config[parts[0]]
            for subkey in parts[1:]:
                value = value[subkey]
            return value
        except KeyError:
            return default


conf = Config.load('.gwfconf.json', defaults=CONFIG_DEFAULTS)
