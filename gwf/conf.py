"""Handling of configuration files."""

import os.path
from ConfigParser import SafeConfigParser

USER_CONFIG_FILE = os.path.expanduser("~/.gwfrc")
LOCAL_CONFIG_FILE = '.gwfrc'

CONFIG_PARSER = SafeConfigParser()
settings = CONFIG_PARSER.read([LOCAL_CONFIG_FILE, USER_CONFIG_FILE])
