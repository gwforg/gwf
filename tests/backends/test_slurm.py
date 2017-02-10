import io
import subprocess
from unittest.mock import ANY, call, mock_open, patch

from gwf import *
from gwf.backends.slurm import (SlurmBackend, _call_sbatch, _call_scancel,
                                _call_squeue, _dump_atomic, _find_exe,
                                _get_live_job_states, _read_json)
from gwf.exceptions import BackendError, NoLogFoundError

from .. import GWFTestCase


class SlurmTestCase(GWFTestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend = SlurmBackend()
        self.config = {}

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
            'gwf.backends.slurm._dump_atomic'
        )


class TestSlurmBackendConfigure(SlurmTestCase):

    def test_does_not_create_log_dir_if_it_already_exists(self):
        self.mock_exists.return_value = True
        self.backend.configure(working_dir='/some/dir',
                               config=self.config)
        self.mock_makedirs.assert_not_called()

    def test_creates_log_dir_if_it_does_not_already_exist(self):
        self.mock_exists.return_value = False
        self.backend.configure(working_dir='/some/dir',
                               config=self.config)
        self.mock_makedirs.assert_called_once_with('/some/dir/.gwf/logs')

    def test_jobdb_is_loaded_from_job_db_file_when_it_exists(self):
        self.mock_read_json.return_value = {"TestTarget": '1000'}
        self.mock_get_live_job_states.return_value = {'1000': 'R'}
        self.backend.configure(working_dir='/some/dir',
                               config=self.config)
        self.mock_read_json.assert_any_call('/some/dir/.gwf/slurm-backend-jobdb.json')
        self.assertDictEqual(self.backend._job_db, {'TestTarget': '1000'})


class TestSlurmBackendGetLiveJobStates(SlurmTestCase):

    def test_live_job_states_are_correctly_parser(self):
        fake_squeue_stdout = '\n'.join([
            "36971043;PD;afterok:36970708,afterok:36970710,afterok:36971042",
            "36971044;R;afterok:36971043",
            "36971045;PD;",
        ])

        fake_squeue_stderr = ''

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

        self.workflow.target('TestTarget1') << outputs('test_output1.txt')
        self.workflow.target('TestTarget2') << outputs('test_output2.txt')
        target3 = self.workflow.target('TestTarget3') <<\
            inputs('test_output1.txt', 'test_output2.txt') <<\
            outputs('test_output3.txt')

        backend = SlurmBackend()
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend,
            config={},
        )

        backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)
        backend.submit(target3, prepared_workflow.dependencies[target3])

        self.mock_call_sbatch.assert_any_call(ANY, ['1000', '2000'])
        self.assertEqual(backend._job_db['TestTarget3'], '3000')
        self.assertEqual(backend._live_job_states['3000'], 'H')

    def test_no_dependency_flag_is_set_if_target_has_no_dependencies(self):
        self.mock_read_json.side_effect = [{}]
        self.mock_get_live_job_states.return_value = {}
        self.mock_call_sbatch.return_value = ('1000', '')

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )

        backend = SlurmBackend()
        backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)
        backend.submit(target, prepared_workflow.dependencies[target])

        self.mock_call_sbatch.assert_any_call(ANY, [])
        self.assertEqual(backend._job_db['TestTarget'], '1000')
        self.assertEqual(backend._live_job_states['1000'], 'H')

    def test_job_script_is_properly_compiled_with_all_supported_options(self):
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )

        backend = SlurmBackend()
        backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

        target = Target(
            name='TestTarget',
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
                'qos': 'somename',
            },
            spec='echo hello world'
        )

        script = backend._compile_script(target)

        self.assertIn('#!/bin/bash', script)
        self.assertIn('#SBATCH --job-name=TestTarget', script)
        self.assertIn('#SBATCH -c 16', script)
        self.assertIn('#SBATCH --mem=16g', script)
        self.assertIn('#SBATCH -t 12:00:00', script)
        self.assertIn('#SBATCH -p normal', script)
        self.assertIn('#SBATCH -A someaccount', script)
        self.assertIn('#SBATCH -C graphics*4', script)
        self.assertIn('#SBATCH --mail-type=BEGIN,END,FAIL', script)
        self.assertIn('#SBATCH --mail-user=test@domain.com', script)
        self.assertIn('#SBATCH --qos=somename', script)
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

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )

        backend = SlurmBackend()
        backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)
        backend.cancel(target)

        self.mock_call_scancel.assert_any_call('1000')

    def test_cancelling_non_running_target_raises_exception(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_live_job_states.return_value = {}

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )

        backend = SlurmBackend()
        backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

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

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

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

        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

        self.assertTrue(self.backend.running(target1))
        self.assertFalse(self.backend.running(target2))
        self.assertFalse(self.backend.running(target3))


class TestSlurmBackendLogs(SlurmTestCase):

    def test_logs_raises_exception_target_has_no_log(self):
        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

        with self.assertRaises(NoLogFoundError):
            self.backend.logs(target)

    def test_logs_returns_log_if_target_has_been_run_once(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
        ]
        self.mock_get_live_job_states.return_value = {
            '1000': 'R'
        }

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = self.backend.logs(target)
            self.assertEqual(stdout.read(), 'this is the log file')
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.stdout', 'r')

    def test_logs_returns_stderr_if_stderr_is_true(self):
        self.mock_read_json.side_effect = [
            {'TestTarget': '1000'},
        ]
        self.mock_get_live_job_states.return_value = {
            '1000': 'R'
        }

        target = self.workflow.target('TestTarget')
        prepared_workflow = PreparedWorkflow(
            targets=self.workflow.targets,
            working_dir=self.workflow.working_dir,
            supported_options=SlurmBackend.supported_options,
            config={},
        )
        self.backend.configure(
            working_dir=prepared_workflow.working_dir, config=self.config)

        m = mock_open()
        m.side_effect = [
            io.StringIO('this is stderr')
        ]

        with patch('builtins.open', m):
            stderr = self.backend.logs(target, stderr=True)
            self.assertEqual(stderr.read(), 'this is stderr')
            m.assert_has_calls([
                call('/some/dir/.gwf/logs/TestTarget.stderr', 'r')
            ])


class TestSlurmBackendClose(SlurmTestCase):

    def test_close_dumps_job_db_if_configure_has_been_called(self):
        self.backend.configure(
            working_dir='/some/dir',
            config=self.config
        )
        self.backend.close()
        self.mock_dump_atomic.assert_any_call(
            ANY,
            '/some/dir/.gwf/slurm-backend-jobdb.json'
        )

    def test_close_does_not_dump_job_db_if_configure_has_not_been_called(self):
        self.backend.close()
        self.mock_dump_atomic.assert_not_called()


class TestFindExe(GWFTestCase):

    def setUp(self):
        self.mock_find_executable = self.create_patch(
            'gwf.backends.slurm.find_executable'
        )

    def test_calls_find_executable_with_name_and_return_result_if_not_none(self):
        self.mock_find_executable.return_value = '/bin/test'

        # Access the decorated function, that is, the actual _find_exe function.
        # If we do not do this, the cache will cache the /bin/test result for
        # 'test' and mess up the remaining tests.
        result = _find_exe.__wrapped__('test')

        self.mock_find_executable.assert_called_once_with('test')
        self.assertEqual(result, '/bin/test')

    def test_raises_exception_if_find_executable_returns_none(self):
        self.mock_find_executable.return_value = None

        with self.assertRaises(BackendError):
            _find_exe.__wrapped__('test')

        self.mock_find_executable.assert_called_once_with('test')


class TestDumpAtomic(GWFTestCase):

    def setUp(self):
        self.mock_open = self.create_patch(
            'builtins.open'
        )
        self.mock_json_dump = self.create_patch(
            'gwf.backends.slurm.json.dump'
        )
        self.mock_os_fsync = self.create_patch(
            'gwf.backends.slurm.os.fsync'
        )
        self.mock_os_rename = self.create_patch(
            'gwf.backends.slurm.os.rename'
        )

    def test_dumps_file_atomically(self):
        mock_fileobj = self.mock_open.return_value.__enter__.return_value
        mock_fileobj.fileno.return_value = 3

        _dump_atomic({'hello': 'world'}, '/some/dir/db.json')

        self.mock_open.assert_called_once_with(
            '/some/dir/db.json.new', 'w'
        )
        self.mock_json_dump.assert_called_once_with(
            {'hello': 'world'}, mock_fileobj
        )
        mock_fileobj.flush.assert_called_once_with()
        self.mock_os_fsync.assert_called_with(3)
        mock_fileobj.close.assert_called_once_with()
        self.mock_os_rename.assert_called_once_with('/some/dir/db.json.new',
                                                    '/some/dir/db.json')


class TestReadJson(GWFTestCase):

    def setUp(self):
        self.mock_open = self.create_patch(
            'builtins.open'
        )
        self.mock_json_load = self.create_patch(
            'gwf.backends.slurm.json.load'
        )

    def test_return_results_of_json_load_if_file_exists(self):
        mock_fileobj = self.mock_open.return_value.__enter__.return_value
        self.mock_json_load.return_value = {'hello': 'world'}

        result = _read_json('test.json')

        self.mock_json_load.assert_called_once_with(mock_fileobj)
        self.assertDictEqual(result, {'hello': 'world'})

    def test_return_empty_dict_if_file_does_not_exist(self):
        self.mock_open.side_effect = OSError

        result = _read_json('test.json')
        self.assertDictEqual(result, {})


class TestCallSqueue(GWFTestCase):

    def setUp(self):
        self.mock_find_exe = self.create_patch(
            'gwf.backends.slurm._find_exe'
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen'
        )

    def test_cmd_uses_executable_found_and_returns_stdout_and_stderr(self):
        self.mock_find_exe.return_value = '/bin/squeue'
        self.mock_popen.return_value.returncode = 0
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        stdout, stderr = _call_squeue()

        self.mock_popen.assert_called_once_with(
            ['/bin/squeue', '--noheader', '--format=%i;%t;%E'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )

        self.assertEqual(stdout, 'this is stdout')
        self.assertEqual(stderr, 'this is stderr')

    def test_raises_exception_with_stderr_if_command_fails(self):
        self.mock_find_exe.return_value = '/bin/squeue'
        self.mock_popen.return_value.returncode = 42
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        with self.assertRaisesRegex(BackendError, 'this is stderr'):
            _call_squeue()


class TestCallScancel(GWFTestCase):

    def setUp(self):
        self.mock_find_exe = self.create_patch(
            'gwf.backends.slurm._find_exe'
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen'
        )

    def test_cmd_uses_executable_found_and_runs_it_correctly(self):
        self.mock_find_exe.return_value = '/bin/scancel'
        self.mock_popen.return_value.returncode = 0
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        _call_scancel('1000')

        self.mock_popen.assert_called_once_with(
            ['/bin/scancel', '-j', '1000'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )

    def test_raises_exception_with_stderr_if_command_fails(self):
        self.mock_find_exe.return_value = '/bin/scancel'
        self.mock_popen.return_value.returncode = 42
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        with self.assertRaisesRegex(BackendError, 'this is stderr'):
            _call_scancel('1000')


class TestCallSbatch(GWFTestCase):

    def setUp(self):
        self.mock_find_exe = self.create_patch(
            'gwf.backends.slurm._find_exe'
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen'
        )

    def test_cmd_uses_executable_found_and_runs_correctly_when_no_dependencies(self):
        self.mock_find_exe.return_value = '/bin/sbatch'
        self.mock_popen.return_value.returncode = 0
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        script = 'this is the script'

        _call_sbatch(script, [])

        self.mock_popen.assert_called_once_with(
            ['/bin/sbatch', '--parsable'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )

        self.mock_popen.return_value.communicate.assert_called_once_with(script)

    def test_cmd_uses_executable_found_and_runs_correctly_with_dependencies(self):
        self.mock_find_exe.return_value = '/bin/sbatch'
        self.mock_popen.return_value.returncode = 0
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        script = 'this is the script'

        _call_sbatch(script, ['1', '2'])

        self.mock_popen.assert_called_once_with(
            ['/bin/sbatch', '--parsable', '--dependency=afterok:1,2'],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
            universal_newlines=True,
        )

        self.mock_popen.return_value.communicate.assert_called_once_with(
            script
        )

    def test_raises_exception_with_stderr_if_command_fails(self):
        self.mock_find_exe.return_value = '/bin/sbatch'
        self.mock_popen.return_value.returncode = 42
        self.mock_popen.return_value.communicate.return_value = (
            'this is stdout', 'this is stderr'
        )

        script = 'this is the script'

        with self.assertRaisesRegex(BackendError, 'this is stderr'):
            _call_sbatch(script, [])
