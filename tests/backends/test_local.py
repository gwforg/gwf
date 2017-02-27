import time
from multiprocessing import Process
from unittest import skip

from gwf.backends.local import Client, Server, State
from gwf.core import Target, Workflow

from .. import GWFTestCase


@skip('Local backend is experimental')
class ServerTestCase(GWFTestCase):

    def setUp(self):
        self.server = Server(num_workers=2)
        self._server = Process(target=self.server.start)
        self._server.start()
        self.addCleanup(self._server.terminate)

        time.sleep(0.5)

        self.client = Client(('', 25000))
        self.addCleanup(self.client.close)

        self.workflow = Workflow(working_dir='/tmp')

        self.target = Target(
            name='TestTarget',
            inputs=[],
            outputs=[],
            options={},
            workflow=self.workflow,
            namespace=None,
            spec='sleep 5',
        )

        self.failing_target = Target(
            name='TestTarget',
            inputs=[],
            outputs=[],
            options={},
            workflow=self.workflow,
            namespace=None,
            spec='sleep 5 && exit 1',
        )


@skip('Local backend is experimental')
class TestServer(ServerTestCase):

    def test_can_get_status_when_no_jobs_have_been_submitted(self):
        result = self.client.status()
        self.assertDictEqual(result, {})

    def test_can_get_status_when_of_submitted_target(self):
        target_id = self.client.submit(self.target)
        time.sleep(0.5)
        result = self.client.status()
        self.assertDictEqual(result, {target_id: State.started})

    def test_submitting_more_targets_than_workers_should_leave_some_targets_pending(self):
        target_id1 = self.client.submit(self.target)
        target_id2 = self.client.submit(self.target)
        target_id3 = self.client.submit(self.target)
        time.sleep(0.5)
        result = self.client.status()
        self.assertDictEqual(result, {
            target_id1: State.started,
            target_id2: State.started,
            target_id3: State.pending,
        })

    def test_target_with_dependency_should_run_after_the_dependency(self):
        target_id1 = self.client.submit(self.target)
        target_id2 = self.client.submit(self.target, deps=[target_id1])
        time.sleep(0.5)
        result = self.client.status()
        self.assertDictEqual(result, {
            target_id1: State.started,
            target_id2: State.pending,
        })
        time.sleep(6.0)
        result = self.client.status()
        self.assertDictEqual(result, {
            target_id1: State.completed,
            target_id2: State.started,
        })
        time.sleep(6.0)
        self.assertDictEqual(result, {
            target_id1: State.completed,
            target_id2: State.completed,
        })

    def test_target_with_two_dependencies_should_run_after_the_dependencies(self):
        pass


@skip('Local backend is experimental')
class TestClient(GWFTestCase):
    pass


@skip('Local backend is experimental')
class TestBackend(GWFTestCase):
    pass
