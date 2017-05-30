import contextlib
import subprocess
from unittest.mock import ANY, patch

from gwf import Target
from gwf.backends import Status, UnknownTargetError, UnknownDependencyError
from gwf.backends.slurm import (SlurmBackend, _call_sbatch, _call_scancel,
                                _call_squeue, _find_exe, init_status_from_queue)
from gwf.exceptions import BackendError
from tests import GWFTestCase


class TestFindExe(GWFTestCase):

    def setUp(self):
        self.mock_find_executable = self.create_patch(
            'gwf.backends.slurm.find_executable', autospec=True
        )

    def test_calls_find_executable_with_name_and_return_result_if_not_none(self):
        self.mock_find_executable.return_value = '/bin/test'

        # Access the decorated function, that is, the actual _find_exe function.
        # If we do not do this, the cache will cache the /bin/test result for
        # 'test' and mess up the remaining tests.
        result = _find_exe.__wrapped__('test')
        self.assertEqual(result, '/bin/test')

    def test_raises_exception_if_find_executable_returns_none(self):
        self.mock_find_executable.return_value = None
        with self.assertRaises(BackendError):
            _find_exe.__wrapped__('test')


class SlurmTestCase(GWFTestCase):

    def setUp(self):
        self.mock_init_status_from_queue = self.create_patch(
            'gwf.backends.slurm.init_status_from_queue',
            autospec=True, return_value={},
        )

        self.mock_call_squeue = self.create_patch(
            'gwf.backends.slurm._call_squeue',
            autospec=True
        )

        self.mock_call_sbatch = self.create_patch(
            'gwf.backends.slurm._call_sbatch',
            autospec=True, return_value=('', '')
        )

        self.mock_call_scancel = self.create_patch(
            'gwf.backends.slurm._call_scancel',
            autospec=True
        )

    @contextlib.contextmanager
    def force_job_id(self, job_id):
        self.mock_call_sbatch.return_value = ('{}\n'.format(job_id), '')
        yield
        self.mock_call_sbatch.return_value = ('', '')


class TestInitStatesFromQueue(SlurmTestCase):

    def test_parse_squeue_output(self):
        fake_squeue_stdout = '\n'.join([
            "36971043;PD;afterok:36970708,afterok:36970710,afterok:36971042",
            "36971044;R;afterok:36971043",
            "36971045;PD;",
        ])

        fake_squeue_stderr = ''

        self.mock_call_squeue.return_value = (
            fake_squeue_stdout, fake_squeue_stderr
        )

        result = init_status_from_queue()

        self.assertDictEqual(result, {
            '36971043': Status.SUBMITTED,
            '36971044': Status.RUNNING,
            '36971045': Status.SUBMITTED,
        })


class TestSlurmBackendSubmit(SlurmTestCase):

    def test_submitting_target_correctly_sets_dependency_flag_for_sbatch(self):
        target1 = Target(
            'TestTarget1',
            inputs=[],
            outputs=['test_output1.txt'],
            options={},
            working_dir='/some/dir',
        )
        target2 = Target(
            'TestTarget2',
            inputs=[],
            outputs=['test_output2.txt'],
            options={},
            working_dir='/some/dir',
        )
        target3 = Target(
            'TestTarget3',
            inputs=['test_output1.txt', 'test_output2.txt'],
            outputs=['test_output3.txt'],
            options={},
            working_dir='/some/dir',
        )

        backend = SlurmBackend()
        with self.force_job_id(1000):
            backend.submit(target1, [])
        with self.force_job_id(2000):
            backend.submit(target2, [])
        with self.force_job_id(3000):
            backend.submit(target3, [target1, target2])

        self.mock_call_sbatch.assert_any_call(ANY, ['1000', '2000'])

    def test_no_dependency_flag_is_set_if_target_has_no_dependencies(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend()
        backend.submit(target, [])

        self.mock_call_sbatch.assert_any_call(ANY, [])

    def test_submit_target_with_unknown_dependency_raises_unknown_dependency_error(self):
        target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend()
        with self.assertRaises(UnknownDependencyError):
            backend.submit(target2, [target1])

    def test_job_script_is_properly_compiled_with_all_supported_options(self):
        backend = SlurmBackend()

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


class TestSlurmBackendCancel(SlurmTestCase):

    def test_cancel_submitted_target(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend()
        with self.force_job_id(1000):
            backend.submit(target, dependencies=[])

        self.assertEqual(backend.status(target), Status.SUBMITTED)
        backend.cancel(target)
        self.assertEqual(backend.status(target), Status.UNKNOWN)

    def test_cancel_unknown_target_raises_exception(self):
        target = Target('TestTarget', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend()
        with self.assertRaises(UnknownTargetError):
            backend.cancel(target)


class TestSlurmBackendSubmitted(SlurmTestCase):

    def test_submitted_should_only_return_true_if_target_is_in_job_db(self):
        target1 = Target('TestTarget1', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        target2 = Target('TestTarget2', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        backend = SlurmBackend()
        with self.force_job_id(1000):
            backend.submit(target1, dependencies=[])

        self.assertEqual(backend.status(target1), Status.SUBMITTED)
        self.assertEqual(backend.status(target2), Status.UNKNOWN)



class TestSlurmBackendClose(SlurmTestCase):

    def test_persist_tracked_targets_on_close(self):
        backend = SlurmBackend()

        with patch.object(backend._tracked, 'persist') as mock_persist:
            backend.close()
            mock_persist.assert_called_once_with()


class TestCallSqueue(GWFTestCase):

    def setUp(self):
        self.mock_find_exe = self.create_patch(
            'gwf.backends.slurm._find_exe', autospec=True
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen', autospec=True
        )

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
            'gwf.backends.slurm._find_exe', autospec=True
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen', autospec=True
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
            'gwf.backends.slurm._find_exe', autospec=True
        )
        self.mock_popen = self.create_patch(
            'gwf.backends.slurm.subprocess.Popen', autospec=True
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
