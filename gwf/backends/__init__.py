from __future__ import absolute_import, print_function

from gwf.backends.slurm import SlurmBackend

AVAILABLE_BACKENDS = {
    'slurm': SlurmBackend,
}
