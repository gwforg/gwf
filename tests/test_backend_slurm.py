import io
import unittest
from unittest.mock import ANY, Mock, call, mock_open, patch

from gwf import PreparedWorkflow, Target, Workflow
from gwf.backends.slurm import SlurmBackend, _get_live_job_states
from gwf.exceptions import BackendError, NoLogFoundError


class SlurmTestCase(unittest.TestCase):

    def create_patch(self, name):
        patcher = patch(name)
        thing = patcher.start()
        self.addCleanup(patcher.stop)
        return thing

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend = SlurmBackend()

        self.mock_get_live_job_states = self.create_patch(
            'gwf.backends.slurm._get_live_job_states'
        )
        self.mock_exists = self.create_patch(
            'gwf.backends.slurm.os.path.exists'
        )
        self.mock_makedirs = self.create_patch(
            'gwf.backends.slurm.os.makedirs'
        )
        self.mock_read_json = self.create_patch(
            'gwf.backends.slurm._read_json'
        )
        self.mock_call_squeue = self.create_patch(
            'gwf.backends.slurm._call_squeue'
        )
        self.mock_call_sbatch = self.create_patch(
            'gwf.backends.slurm._call_sbatch'
        )
        self.mock_call_scancel = self.create_patch(
            'gwf.backends.slurm._call_scancel'
        )
        self.mock_dump_atomic = self.create_patch(
            'gwf.backends.slurm.dump_atomic'
        )


class TestSlurmBackendConfigure(SlurmTestCase):

    def test_does_not_create_log_dir_if_it_already_exists(self):
        self.mock_exists.return_value = True
        self.backend.configure(workflow=self.prepared_workflow)
        self.mock_makedirs.assert_not_called()

    def test_creates_log_dir_if_it_does_not_already_exist(self):
        self.mock_exists.return_value = False
        self.backend.configure(workflow=self.prepared_workflow)
        self.mock_makedirs.assert_called_once_with('/some/dir/.gwf/logs')

    def test_jobdb_is_empty_if_job_db_file_does_not_exist(self):
        self.mock_exists.return_value = True
        self.backend.configure(workflow=self.prepared_workflow)
        self.mock_read_json.assert_any_call(
            '.gwf/slurm-backend-jobdb.json')
        self.assertDictEqual(self.backend._job_db, {})

    def test_jobdb_is_loaded_from_job_db_file_when_it_exists(self):
        self.mock_read_json.return_value = {"TestTarget": '1000'}
        self.mock_get_live_job_states.return_value = {'1000': 'R'}
        self.backend.configure(workflow=self.prepared_workflow)
        self.mock_read_json.assert_any_call('.gwf/slurm-backend-jobdb.json')
        self.assertDictEqual(self.backend._job_db, {'TestTarget': '1000'})


class TestSlurmBackendGetLiveJobStates(SlurmTestCase):

    def test_live_job_states_are_correctly_parser(self):
        fake_squeue_stdout = io.StringIO('\n'.join([
            "36971043;PD;afterok:36970708,afterok:36970710,afterok:36971042",
            "36971044;R;afterok:36971043",
            "36971045;PD;",
        ]))

        fake_squeue_stderr = io.StringIO('')

        self.mock_call_squeue.return_value = (
            fake_squeue_stdout, fake_squeue_stderr
        )

        result = _get_live_job_states()

        self.assertDictEqual(result, {
            '36971043': 'H',
            '36971044': 'R',
            '36971045': 'Q',
        })


class TestSlurmBackendSubmit(SlurmTestCase):

    def test_submitting_target_correctly_sets_dependency_flag_for_sbatch(self):
        self.mock_read_json.side_effect = [
            {'TestTarget1': '1000', 'TestTarget2': '2000'}, {}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'R', '2000': 'H'}
        self.mock_call_sbatch.return_value = ('3000', '')

        target1 = self.workflow.target(
            'TestTarget1',
            outputs=['test_output1.txt'],
        ) << ''
        target2 = self.workflow.target(
            'TestTarget2',
            outputs=['test_output2.txt']
        ) << ''
        target3 = self.workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        ) << ''

        prepared_workflow = PreparedWorkflow(self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)
        backend.submit(target3)

        self.mock_call_sbatch.assert_any_call(ANY, ['1000', '2000'])
        self.assertEqual(backend._job_db['TestTarget3'], '3000')
        self.assertEqual(backend._live_job_states['3000'], 'H')
        self.assertEqual(backend._job_history['TestTarget3'], ['3000'])

    def test_no_dependency_flag_is_set_if_target_has_no_dependencies(self):
        self.mock_read_json.side_effect = [{}, {}]
        self.mock_get_live_job_states.return_value = {}
        self.mock_call_sbatch.return_value = ('1000', '')

        target = self.workflow.target('TestTarget') << ''
        prepared_workflow = PreparedWorkflow(self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)
        backend.submit(target)

        self.mock_call_sbatch.assert_any_call(ANY, [])
        self.assertEqual(backend._job_db['TestTarget'], '1000')
        self.assertEqual(backend._live_job_states['1000'], 'H')
        self.assertEqual(backend._job_history['TestTarget'], ['1000'])

    def test_submitting_target_twice_appends_new_job_id_to_job_history(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'R'}
        self.mock_call_sbatch.return_value = ('2000', '')

        target = self.workflow.target('TestTarget') << ''
        prepared_workflow = PreparedWorkflow(self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)
        backend.submit(target)

        self.assertEqual(backend._job_db['TestTarget'], '2000')
        self.assertEqual(backend._live_job_states['2000'], 'H')
        self.assertEqual(backend._live_job_states['1000'], 'R')
        self.assertEqual(backend._job_history['TestTarget'], ['1000', '2000'])

    def test_job_script_is_properly_compiled_with_all_supported_options(self):
        prepared_workflow = PreparedWorkflow(workflow=self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)

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


class TestSlurmBackendCancel(SlurmTestCase):

    def test_cancelling_a_target_calls_scancel_with_correct_job_id(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'R'}

        target = self.workflow.target('TestTarget') << ''
        prepared_workflow = PreparedWorkflow(self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)
        backend.cancel(target)

        self.mock_call_scancel.assert_any_call('1000')

    def test_cancelling_non_running_target_raises_exception(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)

        backend = SlurmBackend()
        backend.configure(workflow=prepared_workflow)

        with self.assertRaises(BackendError):
            backend.cancel(target)


class TestSlurmBackendSubmitted(SlurmTestCase):

    def test_submitted_should_only_return_true_if_target_is_in_job_db(self):
        self.mock_read_json.side_effect = [
            {'TestTarget1': '1000'},
            {'TestTarget1': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'H'}

        target1 = self.workflow.target('TestTarget1')
        target2 = self.workflow.target('TestTarget2')

        prepared_workflow = PreparedWorkflow(workflow=self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        self.assertTrue(self.backend.submitted(target1))
        self.assertFalse(self.backend.submitted(target2))


class TestSlurmBackendRunning(SlurmTestCase):

    def test_running_should_return_true_if_target_is_in_job_db_and_is_running(self):
        self.mock_read_json.side_effect = [
            {'TestTarget1': '1000', 'TestTarget2': '2000'},
            {'TestTarget1': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'R', '2000': 'H'}

        target1 = self.workflow.target('TestTarget1')
        target2 = self.workflow.target('TestTarget2')
        target3 = self.workflow.target('TestTarget3')

        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        self.assertTrue(self.backend.running(target1))
        self.assertFalse(self.backend.running(target2))
        self.assertFalse(self.backend.running(target3))


class TestSlurmBackendLogs(SlurmTestCase):

    def test_logs_raises_exception_if_rewind_is_0_and_target_has_no_history(self):
        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        with self.assertRaises(NoLogFoundError):
            self.backend.logs(target)

    def test_logs_returns_log_if_rewind_is_0_and_target_has_been_run_once(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {'1000': 'R'}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = self.backend.logs(target)
            self.assertEqual(stdout.read(), 'this is the log file')
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.1000.stdout')

    def test_logs_returns_log_if_rewind_is_1_and_target_has_been_run_twice(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['200', '1000']}
        ]
        self.mock_get_live_job_states.return_value = {}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = self.backend.logs(target, rewind=1)

            self.assertEqual(
                stdout.read(), 'this is the log file'
            )
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.200.stdout'
            )

    def test_logs_raises_exception_if_rewind_is_2_and_target_has_been_run_twice(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['200', '1000']}
        ]
        self.mock_get_live_job_states.return_value = {}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        with self.assertRaises(NoLogFoundError):
            self.backend.logs(target, rewind=2)

    def test_logs_returns_both_stdout_and_stderr_if_stderr_is_true(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(self.workflow)
        self.backend.configure(workflow=prepared_workflow)

        m = mock_open()
        m.side_effect = [
            io.StringIO('this is stdout'),
            io.StringIO('this is stderr')
        ]

        with patch('builtins.open', m):
            stdout, stderr = self.backend.logs(target, stderr=True)
            self.assertEqual(stdout.read(), 'this is stdout')
            self.assertEqual(stderr.read(), 'this is stderr')
            m.assert_has_calls([
                call('/some/dir/.gwf/logs/TestTarget.1000.stdout'),
                call('/some/dir/.gwf/logs/TestTarget.1000.stderr')
            ])


class TestSlurmBackendClose(SlurmTestCase):

    def test_closing_backend_dumps_database_atomically(self):
        self.backend.configure(workflow=self.prepared_workflow)
        self.backend.close()
        self.mock_dump_atomic.assert_any_call(
            {},
            '.gwf/slurm-backend-jobdb.json'
        )
