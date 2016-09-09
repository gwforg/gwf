# -*- coding: utf-8 -*-

"""Handling of configuration files."""

from __future__ import (absolute_import, print_function, division,
                        unicode_literals)

import os.path
from configparser import SafeConfigParser

USER_CONFIG_FILE = os.path.expanduser("~/.gwfrc")
LOCAL_CONFIG_FILE = '.gwfrc'

settings = SafeConfigParser()
settings.read([LOCAL_CONFIG_FILE, USER_CONFIG_FILE])
