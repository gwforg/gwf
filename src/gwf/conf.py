import json
from collections import ChainMap

import attrs

CONFIG_DEFAULTS = {"verbose": "info", "clean_logs": True, "use_spec_hashes": True}

VERBOSITY_LEVELS = ["warning", "debug", "info", "error"]


def try_int(value):
    try:
        return int(value)
    except ValueError:
        return None


def try_true(value):
    if value in ("true", "yes"):
        return True
    return None


def try_false(value):
    if value in ("false", "no"):
        return False
    return None


CONVERTERS = (
    try_int,
    try_true,
    try_false,
    str,
)


def try_conv(value, converters):
    values = (conv(value) for conv in converters)
    return next(filter(lambda x: x is not None, values), None)


@attrs.define
class FileConfig:
    path: str = attrs.field()
    data: ChainMap = attrs.field()

    def get(self, key, default=None):
        try:
            return self.data[key]
        except KeyError:
            return default

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        print(try_conv(value.strip(), CONVERTERS))
        self.data[key] = try_conv(value, CONVERTERS)

    def __delitem__(self, key):
        if key in self.data:
            del self.data[key]

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def items(self):
        return self.data.items()

    def get_namespace(self, ns):
        res = {}
        for k, v in self.items():
            if k.startswith(ns):
                new_key = k[len(ns) + 1 :]
                res[new_key] = v
        return res

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
