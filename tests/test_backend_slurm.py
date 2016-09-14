import unittest
from unittest.mock import mock_open, patch

from gwf import PreparedWorkflow, Target, Workflow
from gwf.backends.slurm import SlurmBackend


@patch('gwf.backends.slurm.subprocess.Popen')
@patch('gwf.backends.slurm._find_exe', side_effect=['squeue', 'sbatch', 'scancel'])
class TestSlurmBackend(unittest.TestCase):

    def test_job_script_is_properly_compile_with_all_options(self, mock_find_exe, mock_popen):
        backend = SlurmBackend(PreparedWorkflow())

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

        script = backend._compile_script(target)

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

    def test_jobdb_is_empty_if_job_db_file_does_not_exist(self, mock_find_exe, mock_popen):
        workflow = Workflow()
        prepared_workflow = PreparedWorkflow(workflow)

        with patch('builtins.open', side_effect=FileNotFoundError) as mock_open_:
            backend = SlurmBackend(prepared_workflow)
            backend.configure()

            mock_open_.assert_called_once_with('.gwf/slurm-backend-jobdb.json')
            self.assertDictEqual(backend._job_db, {})

    def test_jobdb_is_loaded_from_job_db_file_when_it_exists(self, mock_find_exe, mock_popen):
        workflow = Workflow()
        prepared_workflow = PreparedWorkflow(workflow)

        with patch('builtins.open', mock_open(read_data='{"TestTarget": 1000}')) as mock_open_:
            backend = SlurmBackend(prepared_workflow)
            with patch.object(backend, '_live_job_states') as mock_live_job_states:
                mock_live_job_states.return_value = {1000: 'R'}

                backend.configure()
                mock_open_.assert_called_once_with(
                    '.gwf/slurm-backend-jobdb.json'
                )
                self.assertDictEqual(backend._job_db, {'TestTarget': 'R'})
