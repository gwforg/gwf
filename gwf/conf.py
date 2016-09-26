"""Handling of configuration files."""

import json
import os.path
from collections import ChainMap
from collections.abc import MutableMapping

USER_CONFIG_FILE = os.path.expanduser("~/.gwfrc")
LOCAL_CONFIG_FILE = '.gwfrc'


class Config(MutableMapping):

    def __init__(self, path):
        super().__init__()
        self.path = path
        try:
            with open(self.path) as fileobj:
                self.dct = json.load(fileobj)
        except IOError:
            self.dct = {}

    def _find_option(self, dct, name):
        *segments, last_segment = name.split('.')
        for segment in segments:
            if segment not in dct:
                dct[segment] = {}
            dct = dct[segment]
        return dct, last_segment

    def __setitem__(self, name, value):
        dct, last_segment = self._find_option(self.dct, name)
        dct[last_segment] = value

    def __delitem__(self, name):
        dct, last_segment = self._find_option(self.dct, name)
        del dct[last_segment]

    def __getitem__(self, name):
        dct, last_segment = self._find_option(self.dct, name)
        return dct[last_segment]

    def __iter__(self):
        return iter(self.dct)

    def __len__(self):
        return len(self.dct)

    def __repr__(self):
        return repr(self.dct)

    def __str__(self):
        return str(self.dct)

    def dump(self):
        with open(self.path, 'w') as fileobj:
            json.dump(self.dct, fileobj)


local_settings = Config(LOCAL_CONFIG_FILE)
user_settings = Config(USER_CONFIG_FILE)

defaults = {
    'backend': 'testing'
}

settings = ChainMap(local_settings, user_settings, defaults)
