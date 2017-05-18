import json
from collections import ChainMap


class Config:

    def __init__(self, path, defaults=None, override=None):
        self.path = path
        configs = []
        try:
            with open(path) as config_file:
                config = json.load(config_file)
        except FileNotFoundError:
            config = {}

        if defaults is None:
            defaults = {}
        configs.append(defaults)
        configs.append(config)
        if override is not None:
            configs.append(override)
        self._config = ChainMap(*configs[::-1])

    def get_config(self):
        return dict(self._config)

    def get_file_config(self):
        return Config(path=self.path, defaults=self._config.maps[1])

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

    def get(self, key, default=None):
        parts = key.split('.')
        try:
            value = self._config[parts[0]]
            for subkey in parts[1:]:
                value = value[subkey]
            return value
        except KeyError:
            return default
