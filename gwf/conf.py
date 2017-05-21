import json
from collections import ChainMap


CONFIG_DEFAULTS = {
    'verbose': 'info',
    'backend': 'local',
    'no_color': False,
}


def tokenize_key(key):
    return key.split('.')


class Config:

    def __init__(self, initial):
        self._config = initial

    def __getitem__(self, key):
        tokens = tokenize_key(key)
        value = self._config[tokens[0]]
        for token in tokens[1:]:
            if not isinstance(value, dict):
                raise KeyError(key)
            value = value[token]
        return value

    def __setitem__(self, key, new_value):
        tokens = tokenize_key(key)
        value = self._config
        for token in tokens[:-1]:
            if token not in value:
                value[token] = {}
            value = value[token]
        value[tokens[-1]] = new_value

    def __delitem__(self, key):
        tokens = tokenize_key(key)
        value = self._config
        for token in tokens[:-1]:
            value = value[token]
        del value[tokens[-1]]

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def __len__(self):
        return len(self._config)

    def iterkeys(self):
        def deepiter(dct, parents):
            for key, value in dct.items():
                new_parents = parents + [key]
                key_path = '.'.join(new_parents)
                if isinstance(value, dict):
                    yield from deepiter(value, parents=new_parents)
                else:
                    yield key_path
        yield from deepiter(self._config, parents=[])

    def __iter__(self):
        return self.iterkeys()

    def __repr__(self):
        return repr(self._config)


class FileConfig(Config):

    def __init__(self, path, initial):
        super().__init__(initial)
        self.path = path

    def dump(self):
        """Dump the configuration to disk."""
        with open(self.path, 'w') as config_file:
            json.dump(dict(self._config), config_file, indent=4, sort_keys=True)

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
config = ChainMap(file_config, Config(CONFIG_DEFAULTS))
