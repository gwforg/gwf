import unittest
from unittest.mock import patch

from gwf import PreparedWorkflow, Target, Workflow
from gwf.backends.slurm import SlurmBackend, _compile_script


@patch('gwf.backends.slurm.subprocess.Popen')
@patch('gwf.backends.slurm._find_exe', side_effect=['squeue', 'sbatch', 'scancel'])
class TestSlurmBackend(unittest.TestCase):

    def test_job_script_is_properly_compile_with_all_options(self, mock_find_exe, mock_popen):
        target = Target(
            name='TestTarget',
            inputs=[],
            outputs=[],
            working_dir='/some/dir',
            options={
                'cores': 16,
                'memory': '16g',
                'walltime': '12:00:00',
                'queue': 'normal',
                'account': 'someaccount',
                'constraint': 'graphics*4',
                'mail_type': 'BEGIN,END,FAIL',
                'mail_user': 'test@domain.com'
            },
            spec='echo hello world'
        )

        script = _compile_script(target)

        self.assertIn('#!/bin/bash', script)
        self.assertIn('#SBATCH -c 16', script)
        self.assertIn('#SBATCH --mem=16g', script)
        self.assertIn('#SBATCH -t 12:00:00', script)
        self.assertIn('#SBATCH -p normal', script)
        self.assertIn('#SBATCH -A someaccount', script)
        self.assertIn('#SBATCH -C graphics*4', script)
        self.assertIn('#SBATCH --mail-type=BEGIN,END,FAIL', script)
        self.assertIn('#SBATCH --mail-user=test@domain.com', script)
        self.assertIn('cd /some/dir', script)
        self.assertIn('export GWF_JOBID=$SLURM_JOBID', script)
        self.assertIn('echo hello world', script)
