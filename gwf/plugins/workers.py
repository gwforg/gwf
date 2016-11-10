import multiprocessing

from . import Plugin
from ..backends.local import Server


class WorkersPlugin(Plugin):

    def setup_argument_parser(self, parser, subparsers):
        subparser = self.setup_subparser(
            subparsers,
            'workers',
            'Start worker.s',
            self.on_start_workers,
        )

        subparser.add_argument(
            '-p',
            '--port',
            default=12345,
            type=int,
            help='Port where worker manager accepts clients.',
        )

        subparser.add_argument(
            '-n',
            '--num-workers',
            default=multiprocessing.cpu_count(),
            type=int,
            help='Number of workers to spawn.',
        )

    def on_start_workers(self):
        port = self.config['port']
        num_workers = self.config['num_workers']
        server = Server(port=port, num_workers=num_workers)
        server.start()
