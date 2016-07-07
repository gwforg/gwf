from __future__ import absolute_import
from __future__ import print_function

from gwf.backends.local import LocalBackend
from gwf.backends.slurm import SlurmBackend
from gwf.backends.torque import TorqueBackend


# This will be set in the gwf script and refer to the grid backend used.
BACKEND = None

AVAILABLE_BACKENDS = {
    'torque': TorqueBackend,
    'slurm': SlurmBackend,
    'local': LocalBackend
}
