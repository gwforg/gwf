import unittest
from gwf.cli.arg_parsing import SubCommand, ArgumentDispatching
from gwf.cli.status import StatusCommand
from gwf.exceptions import GWFError, TargetDoesNotExistsError

class TestHandler(SubCommand):
    def __init__(self):
        self.was_called = False

    def set_arguments(self, parser):
        parser.add_argument("--test", action = "store_true")


    def handle(self, arguments):
        self.was_called = True
        self.arguments = arguments


class TestCLI(unittest.TestCase):

    def test_basic_argument_parsing(self):
        dispatcher = ArgumentDispatching()
        test_handler = TestHandler()
        arguments = ["test"]
        dispatcher.install_subcommand("test", "", test_handler)
        dispatcher.dispatch(arguments)
        self.assertTrue(test_handler.was_called)

    def test_subcommand_arguments(self):
        dispatcher = ArgumentDispatching()
        test_handler = TestHandler()
        dispatcher.install_subcommand("test", "", test_handler)

        arguments = ["test"]
        dispatcher.dispatch(arguments)
        self.assertFalse(test_handler.arguments.test)

        arguments = ["test", "--test"]
        dispatcher.dispatch(arguments)
        self.assertTrue(test_handler.arguments.test)

    def test_access_to_workflow(self):
        dispatcher = ArgumentDispatching()
        test_handler = TestHandler()
        dispatcher.install_subcommand("test", "", test_handler)
        dispatcher.dispatch(["test"])
        with self.assertRaises(GWFError):
            test_handler.workflow

        test_handler = TestHandler()
        dispatcher.install_subcommand("test", "", test_handler)
        dispatcher.dispatch(["-w", "examples/minimal-workflow/workflow.py", "test"])
        workflow = test_handler.workflow


class MockStatusCommand(StatusCommand):

    def __init__(self):
        self.targets = []
        self.verbose_called = False
        self.progress_called = False
        super(MockStatusCommand, self).__init__()

    def print_verbose(self, target_names):
        self.verbose_called = True
        self.targets = target_names

    def print_progress(self, target_names):
        self.progress_called = True
        self.targets = target_names

class TestStatusSubcommand(unittest.TestCase):

    def setUp(self):
        self.dispatcher = ArgumentDispatching()
        self.handler = MockStatusCommand()
        self.dispatcher.install_subcommand("status", "", self.handler)


    def test_status_verbose(self):
        self.dispatcher.dispatch(["-w", "examples/minimal-workflow/workflow.py", "status", "--verbose"])
        self.assertTrue(self.handler.verbose_called)
        self.assertFalse(self.handler.progress_called)
        self.assertListEqual(self.handler.targets, ['All'])

    def test_status_progress(self):
        self.dispatcher.dispatch(["-w", "examples/minimal-workflow/workflow.py", "status"])
        self.assertFalse(self.handler.verbose_called)
        self.assertTrue(self.handler.progress_called)
        self.assertListEqual(self.handler.targets, ['All'])

    def test_status_with_named_targets(self):
        self.dispatcher.dispatch(["-w", "examples/minimal-workflow/workflow.py", "status", "World", "Universe"])
        self.assertFalse(self.handler.verbose_called)
        self.assertTrue(self.handler.progress_called)
        self.assertListEqual(self.handler.targets, ["World", "Universe"])

    def test_exception_with_nonexisting_target(self):
        with self.assertRaises(TargetDoesNotExistsError):
            self.dispatcher.dispatch(["-w", "examples/minimal-workflow/workflow.py", "status", "nonexisting_target"])
