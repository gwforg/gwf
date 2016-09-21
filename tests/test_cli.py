import unittest
from gwf.cli.arg_parsing import SubCommand, ArgumentDispatching
from gwf.exceptions import GWFError

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
