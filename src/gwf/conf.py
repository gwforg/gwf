import json
from collections import ChainMap

import attrs

CONFIG_DEFAULTS = {"verbose": "info", "clean_logs": True, "use_spec_hashes": True}


@attrs.define
class FileConfig:
    path: str = attrs.field()
    data: ChainMap = attrs.field()
    validators: dict = attrs.field(init=False, factory=dict)

    def validator(self, key):
        """Register a configuration key validator function."""

        def _inner(func):
            self.validators[key] = func

        return _inner

    def _validate_value(self, key, value):
        if key in self.validators:
            self.validators[key](value)

    def get(self, key, default=None):
        try:
            return self.data[key]
        except KeyError:
            return default

    def __getitem__(self, key):
        value = self.data[key]
        self._validate_value(key, value)
        return value

    def __setitem__(self, key, value):
        self._validate_value(key, value)
        self.data[key] = value

    def __delitem__(self, key):
        if key in self.data:
            del self.data[key]

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def dump(self):
        """Dump the configuration to disk."""
        with open(str(self.path), "w+") as config_file:
            json.dump(dict(self.data.maps[0]), config_file, indent=4, sort_keys=True)

    @classmethod
    def load(cls, path):
        """Load configuration from a file.

        Reads configuration from `file` and returns a :class:`Config` instance
        with the configuration. The `defaults` will be merged into the
        configuration.

        :param path str: Path to the configuration file.
        """
        try:
            with open(str(path)) as config_file:
                data = json.load(config_file)
        except FileNotFoundError:
            data = {}
        return cls(path=path, data=ChainMap(data, CONFIG_DEFAULTS))


def config_from_path(path):
    return FileConfig.load(path)


config = config_from_path(".gwfconf.json")
