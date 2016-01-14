""" Handling of configuration files. """

import os.path
import ConfigParser
from distutils.spawn import find_executable as has_exe


USER_CONFIG_FILE = os.path.expanduser("~/.gwfrc")
LOCAL_CONFIG_FILE = '.gwfrc'


def set_defaults(config):
	"""Set global defaults in the 'gwf' section."""
	config.add_section('gwf')
        if has_exe("squeue") and has_exe("sbatch"):
            config.set('gwf', 'backend', 'slurm')
        else:
            config.set('gwf', 'backend', 'torque')


def read_configurations():
	config = ConfigParser.SafeConfigParser()
	set_defaults(config)
	config.read([USER_CONFIG_FILE, LOCAL_CONFIG_FILE])
	return config


def read_local_configurations():
	"""Read the local configurations only and does not set defaults."""
	config = ConfigParser.SafeConfigParser()
	config.read([LOCAL_CONFIG_FILE])
	return config


def write_local_configurations(config):
	with open(LOCAL_CONFIG_FILE, 'w') as config_file:
		config.write(config_file)


def read_user_configurations():
	"""Read the global/user configurations only and does not set defaults."""
	config = ConfigParser.SafeConfigParser()
	config.read([USER_CONFIG_FILE])
	return config


def write_user_configurations(config):
	with open(USER_CONFIG_FILE, 'w') as config_file:
		config.write(config_file)
