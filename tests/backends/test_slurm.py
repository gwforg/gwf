import io
import subprocess
from unittest.mock import ANY, call, mock_open, patch

from gwf import Graph, Target, Workflow
from gwf.backends.slurm import (SlurmBackend, _call_sbatch, _call_scancel,
                                _call_squeue, _find_exe, _get_status)
from gwf.exceptions import BackendError, NoLogFoundError
from tests import GWFTestCase


class SlurmTestCase(GWFTestCase):

    def setUp(self):
        self.workflow = Workflow(working_dir='/some/dir')
        self.graph = Graph(targets=self.workflow.targets)

        self.mock_get_status = self.create_patch(
            'gwf.backends.slurm._get_status'
        )
        self.mock_exists = self.create_patch(
            'gwf.backends.slurm.os.path.exists'
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
        self.mock_persistabledict = self.create_patch(
            'gwf.backends.slurm.PersistableDict'
        )


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

        result = _get_status()

        self.assertDictEqual(result, {
            '36971043': 'H',
            '36971044': 'R',
            '36971045': 'Q',
        })


class TestSlurmBackendSubmit(SlurmTestCase):

    def test_submitting_target_correctly_sets_dependency_flag_for_sbatch(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget1': '1000', 'TestTarget2': '2000'}, {}
        ]
        self.mock_get_status.return_value = {'1000': 'R', '2000': 'H'}
        self.mock_call_sbatch.return_value = ('3000', '')

        self.workflow.target(
            'TestTarget1',
            inputs=[],
            outputs=['test_output1.txt'],
        )
        self.workflow.target(
            'TestTarget2',
            inputs=[],
            outputs=['test_output2.txt']
        )
        target3 = self.workflow.target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt']
        )

        backend = SlurmBackend(working_dir='/some/dir')
        graph = Graph(targets=self.workflow.targets,)
        backend.submit(target3, graph.dependencies[target3])

        self.mock_call_sbatch.assert_any_call(ANY, ['1000', '2000'])
        self.assertEqual(backend._tracked['TestTarget3'], '3000')
        self.assertEqual(backend._status['3000'], 'H')

    def test_no_dependency_flag_is_set_if_target_has_no_dependencies(self):
        self.mock_persistabledict.side_effect = [{}]
        self.mock_get_status.return_value = {}
        self.mock_call_sbatch.return_value = ('1000', '')

        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        graph = Graph(targets={'TestTarget': target})

        backend = SlurmBackend(working_dir='/some/dir')
        backend.submit(target, graph.dependencies[target])

        self.mock_call_sbatch.assert_any_call(ANY, [])
        self.assertEqual(backend._tracked['TestTarget'], '1000')
        self.assertEqual(backend._status['1000'], 'H')

    def test_job_script_is_properly_compiled_with_all_supported_options(self):
        backend = SlurmBackend(working_dir='/some/dor')

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
        self.mock_persistabledict.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_status.return_value = {'1000': 'R'}

        target = self.workflow.target('TestTarget', inputs=[], outputs=[], working_dir='/some/dir')
        backend = SlurmBackend(working_dir='/some/dir')
        backend.cancel(target)

        self.mock_call_scancel.assert_any_call('1000')

    def test_cancelling_non_running_target_raises_exception(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget': '1000'},
            {'TestTarget': ['1000']}
        ]
        self.mock_get_status.return_value = {}

        target = self.workflow.target('TestTarget', inputs=[], outputs=[], working_dir='/some/dir')
        graph = Graph(targets=self.workflow.targets)
        backend = SlurmBackend(working_dir='/some/dir')
        with self.assertRaises(BackendError):
            backend.cancel(target)


class TestSlurmBackendSubmitted(SlurmTestCase):

    def test_submitted_should_only_return_true_if_target_is_in_job_db(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget1': '1000'},
            {'TestTarget1': ['1000']}
        ]
        self.mock_get_status.return_value = {'1000': 'H'}

        target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend(working_dir='/some/dir')

        self.assertTrue(backend.submitted(target1))
        self.assertFalse(backend.submitted(target2))


class TestSlurmBackendRunning(SlurmTestCase):

    def test_running_should_return_true_if_target_is_in_job_db_and_is_running(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget1': '1000', 'TestTarget2': '2000'},
            {'TestTarget1': ['1000']}
        ]
        self.mock_get_status.return_value = {'1000': 'R', '2000': 'H'}

        target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        target3 = Target('TestTarget3', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend(working_dir='/some/dir')

        self.assertTrue(backend.running(target1))
        self.assertFalse(backend.running(target2))
        self.assertFalse(backend.running(target3))


class TestSlurmBackendLogs(SlurmTestCase):

    def test_logs_raises_exception_target_has_no_log(self):
        target = self.workflow.target('TestTarget', inputs=[], outputs=[], working_dir='/some/dir')
        graph = Graph(targets=self.workflow.targets)
        backend = SlurmBackend(working_dir='/some/dir')
        with self.assertRaises(NoLogFoundError):
            backend.logs(target)

    def test_logs_returns_log_if_target_has_been_run_once(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget': '1000'},
        ]
        self.mock_get_status.return_value = {
            '1000': 'R'
        }

        target = self.workflow.target('TestTarget', inputs=[], outputs=[])
        graph = Graph(targets=self.workflow.targets)
        backend = SlurmBackend(working_dir='/some/dir')

        m = mock_open(read_data='this is the log file')
        with patch('builtins.open', m):
            stdout = backend.logs(target)
            self.assertEqual(stdout.read(), 'this is the log file')
            m.assert_called_once_with(
                '/some/dir/.gwf/logs/TestTarget.stdout', 'r')

    def test_logs_returns_stderr_if_stderr_is_true(self):
        self.mock_persistabledict.side_effect = [
            {'TestTarget': '1000'},
        ]
        self.mock_get_status.return_value = {
            '1000': 'R'
        }

        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        m = mock_open()
        m.side_effect = [
            io.StringIO('this is stderr')
        ]

        backend = SlurmBackend(working_dir=self.workflow.working_dir)
        with patch('builtins.open', m):
            stderr = backend.logs(target, stderr=True)
            self.assertEqual(stderr.read(), 'this is stderr')
            m.assert_has_calls([
                call('/some/dir/.gwf/logs/TestTarget.stderr', 'r')
            ])


class TestSlurmBackendClose(SlurmTestCase):

    def test_persist_tracked_on_close(self):
        backend = SlurmBackend(working_dir=self.workflow.working_dir)
        backend.close()
        self.mock_persistabledict.return_value.persist.assert_called_once_with()


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
            ['/bin/sbatch', '--parsable', '--dependency=afterok:1:2'],
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
