import json
from collections import ChainMap


CONFIG_DEFAULTS = {"verbose": "info", "backend": "local", "no_color": False}


class FileConfig:
    def __init__(self, path, data=None):
        self.path = path

        self._data = ChainMap({}, data)
        self._validators = {}

    def validator(self, key):
        """Register a configuration key validator function."""

        def _inner(func):
            self._validators[key] = func

        return _inner

    def _validate_value(self, key, value):
        if key in self._validators:
            self._validators[key](value)

    def get(self, key, default=None):
        value = self._data.get(key, default)
        self._validate_value(key, value)
        return value

    def __getitem__(self, key):
        value = self._data[key]
        self._validate_value(key, value)
        return value

    def __setitem__(self, key, value):
        self._validate_value(key, value)
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def dump(self):
        """Dump the configuration to disk."""
        with open(self.path, "w") as config_file:
            json.dump(dict(self._data), config_file, indent=4, sort_keys=True)

    @classmethod
    def load(cls, path, data=None):
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

        if data is None:
            data = {}
        else:
            data = dict(data)
        data.update(file_config)
        return cls(path=path, data=data)


config = FileConfig.load(".gwfconf.json", data=CONFIG_DEFAULTS)
