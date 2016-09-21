import io
import subprocess
import unittest
from unittest.mock import Mock, call, mock_open, patch

from gwf import PreparedWorkflow, Target, Workflow
from gwf.backends.slurm import SlurmBackend
from gwf.exceptions import BackendError, NoLogFoundError


@patch('gwf.backends.slurm.subprocess.Popen')
@patch('gwf.backends.slurm._find_exe', side_effect=['squeue', 'sbatch', 'scancel'])
@patch('gwf.backends.slurm.os.makedirs')
class TestSlurmBackend(unittest.TestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.prepared_workflow = PreparedWorkflow(self.workflow)

    @patch('gwf.backends.slurm.os.path.exists', return_value=True)
    def test_does_not_create_log_dir_if_it_already_exists(self, mock_exists, mock_makedirs, mock_find_exe, mock_popen):
        backend = SlurmBackend(self.prepared_workflow)
        backend.configure()
        mock_makedirs.assert_not_called()

    @patch('gwf.backends.slurm.os.path.exists', return_value=False)
    def test_creates_log_dir_if_it_does_not_already_exist(self, mock_exists, mock_makedirs, mock_find_exe, mock_popen):
        backend = SlurmBackend(self.prepared_workflow)
        backend.configure()
        mock_makedirs.assert_called_once_with('/some/dir/.gwf/logs')

    def test_job_script_is_properly_compiled_with_all_supported_options(self, mock_makedirs, mock_find_exe, mock_popen):
        backend = SlurmBackend(self.prepared_workflow)
        backend.configure()

        target = Target(
            name='TestTarget',
            inputs=[],
            outputs=[],
            workflow=self.workflow,
            options={
                'cores': 16,
                'memory': '16g',
                'walltime': '12:00:00',
                'queue': 'normal',
                'account': 'someaccount',
                'constraint': 'graphics*4',
                'mail_type': 'BEGIN,END,FAIL',
                'mail_user': 'test@domain.com',
                'unsupported_option': 'unsupported value'
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

        # Now check that the unsupported option never made it into the script.
        self.assertNotIn('unsupported value', script)

    def test_jobdb_is_empty_if_job_db_file_does_not_exist(self, mock_makedirs, mock_find_exe, mock_popen):
        with patch('builtins.open', side_effect=FileNotFoundError) as mock_open_:
            backend = SlurmBackend(self.prepared_workflow)
            backend.configure()

            mock_open_.assert_any_call('.gwf/slurm-backend-jobdb.json')
            self.assertDictEqual(backend._job_db, {})

    def test_jobdb_is_loaded_from_job_db_file_when_it_exists(self, mock_makedirs, mock_find_exe, mock_popen):
        with patch('builtins.open', mock_open(read_data='{"TestTarget": 1000}')) as mock_open_:
            backend = SlurmBackend(self.prepared_workflow)
            with patch.object(
                    backend,
                    '_get_live_job_states',
                    return_value={1000: 'R'}):

                backend.configure()
                mock_open_.assert_any_call(
                    '.gwf/slurm-backend-jobdb.json'
                )
                self.assertDictEqual(backend._job_db, {'TestTarget': 'R'})

    def test_live_job_states_are_correctly_parser(self, mock_makedirs, mock_find_exe, mock_popen):
        backend = SlurmBackend(self.prepared_workflow)
        backend.configure()

        fake_squeue_output = io.StringIO('\n'.join([
            "36971043;PD;afterok:36970708,afterok:36970710,afterok:36971042",
            "36971044;R;afterok:36971043",
            "36971045;PD;",
        ]))

        mock_popen_instance = Mock(spec=subprocess.Popen)
        mock_popen_instance.stdout = fake_squeue_output
        mock_popen.return_value = mock_popen_instance

        result = backend._get_live_job_states()

        self.assertDictEqual(result, {
            '36971043': 'H',
            '36971044': 'R',
            '36971045': 'Q',
        })

    @patch('gwf.backends.slurm.dump_atomic')
    def test_closing_backend_dumps_database_atomically(self, mock_dump_atomic, mock_makedirs, mock_find_exe, mock_popen):
        backend = SlurmBackend(self.prepared_workflow)
        backend.configure()
        backend.close()

        mock_dump_atomic.assert_any_call({}, '.gwf/slurm-backend-jobdb.json')

    def test_submitted_should_return_true_if_target_is_in_job_db(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target1 = workflow.target('TestTarget1')
        target2 = workflow.target('TestTarget2')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend._job_db = {'TestTarget1': '1000'}

        self.assertTrue(backend.submitted(target1))
        self.assertFalse(backend.submitted(target2))

    def test_running_should_return_true_if_target_is_in_job_db_and_is_running(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target1 = workflow.target('TestTarget1')
        target2 = workflow.target('TestTarget2')
        target3 = workflow.target('TestTarget3')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend._job_db = {'TestTarget1': '1000', 'TestTarget2': '2000'}
        backend._live_job_states = {'1000': 'R', '2000': 'H'}

        self.assertTrue(backend.running(target1))
        self.assertFalse(backend.running(target2))
        self.assertFalse(backend.running(target3))

    def test_submitting_target_correctly_sets_dependency_flag_for_sbatch(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target1 = workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt'],
        ) << ''
        target2 = workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        ) << ''
        target3 = workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        ) << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        backend._job_db = {'TestTarget1': '1000', 'TestTarget2': '2000'}
        backend._live_job_states = {}

        mock_popen_instance = Mock()
        mock_popen_instance.returncode = 0
        mock_popen_instance.communicate.return_value = ('3000\n', '')
        mock_popen.return_value = mock_popen_instance

        backend.submit(target3)

        mock_popen.assert_any_call([
            'sbatch',
            '--parsable',
            '--dependency=afterok:1000,2000',
        ])

        self.assertEqual(backend._job_db['TestTarget3'], '3000')
        self.assertEqual(backend._live_job_states['3000'], 'H')

    def test_no_dependency_flag_is_set_if_target_has_no_dependencies(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget') << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        mock_popen_instance = Mock()
        mock_popen_instance.returncode = 0
        mock_popen_instance.communicate.return_value = ('3000\n', '')
        mock_popen.return_value = mock_popen_instance

        backend.submit(target)

        mock_popen.assert_any_call([
            'sbatch',
            '--parsable',
        ])

    def test_submitting_target_adds_job_id_to_job_history(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget') << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        mock_popen_instance = Mock()
        mock_popen_instance.returncode = 0
        mock_popen_instance.communicate.return_value = ('1000\n', '')
        mock_popen.return_value = mock_popen_instance

        backend.submit(target)
        self.assertIn('TestTarget', backend._job_history)
        self.assertIn('1000', backend._job_history['TestTarget'])

    def test_submitting_target_twice_appends_new_job_id_to_job_history(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget') << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        backend._job_history['TestTarget'] = ['1000']

        mock_popen_instance = Mock()
        mock_popen_instance.returncode = 0
        mock_popen_instance.communicate.return_value = ('2000\n', '')
        mock_popen.return_value = mock_popen_instance

        backend.submit(target)
        self.assertIn('TestTarget', backend._job_history)
        self.assertEqual(backend._job_history['TestTarget'], ['1000', '2000'])

    def test_submitting_target_raises_exception_if_sbatch_cannot_be_called(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target1 = workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt']
        ) << ''
        target2 = workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        ) << ''
        target3 = workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        ) << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()
        backend._job_db = {'TestTarget1': '1000', 'TestTarget2': '2000'}
        backend._live_job_states = {}

        mock_popen_instance = Mock()
        mock_popen_instance.returncode = 1
        mock_popen_instance.communicate.return_value = ('', '')
        mock_popen.return_value = mock_popen_instance

        with self.assertRaises(BackendError):
            backend.submit(target3)

    def test_cancelling_a_target_calls_scancel_with_correct_job_id(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget') << ''

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        backend._job_db = {'TestTarget': '1000'}
        backend._live_job_states = {'1000': 'R'}

        backend.cancel(target)

        mock_popen.assert_any_call(['scancel', '-j', '1000'])

    def test_cancelling_non_running_target_raises_exception(self,  mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        backend._job_db = {'TestTarget': '1000'}
        backend._live_job_states = {}

        with self.assertRaises(BackendError):
            backend.cancel(target)

    def test_logs_raises_exception_if_rewind_is_0_and_target_has_no_history(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow()
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()

        with self.assertRaises(NoLogFoundError):
            backend.logs(target)

    def test_logs_returns_log_if_rewind_is_0_and_target_has_been_run_once(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()
        backend._job_history[target.name] = [1000]

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = backend.logs(target)
            self.assertEqual(stdout.read(), 'this is the log file')
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.1000.stdout')

    def test_logs_returns_log_if_rewind_is_1_and_target_has_been_run_twice(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()
        backend._job_history[target.name] = [200, 1000]

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = backend.logs(target, rewind=1)
            self.assertEqual(stdout.read(), 'this is the log file')
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.200.stdout')

    def test_logs_raises_exception_if_rewind_is_2_and_target_has_been_run_twice(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()
        backend._job_history[target.name] = [200, 1000]

        with self.assertRaises(NoLogFoundError):
            backend.logs(target, rewind=2)

    def test_logs_returns_both_stdout_and_stderr_if_stderr_is_true(self, mock_makedirs, mock_find_exe, mock_popen):
        workflow = Workflow(working_dir='/some/dir')
        target = workflow.target('TestTarget')

        prepared_workflow = PreparedWorkflow(workflow)

        backend = SlurmBackend(prepared_workflow)
        backend.configure()
        backend._job_history[target.name] = [1000]

        m = mock_open()
        m.side_effect = [io.StringIO('this is stdout'),
                         io.StringIO('this is stderr')]
        with patch('builtins.open', m):
            stdout, stderr = backend.logs(target, stderr=True)
            self.assertEqual(stdout.read(), 'this is stdout')
            self.assertEqual(stderr.read(), 'this is stderr')
            m.assert_has_calls([
                call('/some/dir/.gwf/logs/TestTarget.1000.stdout'),
                call('/some/dir/.gwf/logs/TestTarget.1000.stderr')
            ])
