"""Handling of configuration files."""

import os.path
from configparser import SafeConfigParser

USER_CONFIG_FILE = os.path.expanduser("~/.gwfrc")
LOCAL_CONFIG_FILE = '.gwfrc'

settings = SafeConfigParser()
settings.read([LOCAL_CONFIG_FILE, USER_CONFIG_FILE])
