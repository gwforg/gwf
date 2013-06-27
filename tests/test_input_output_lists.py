import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class ListTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'input-output-lists.gwf'))

    def test_lists(self):
        self.assertItemsEqual(['infiles','outfiles'], self.workflow.lists)
        
    def test_target(self):
        self.assertItemsEqual(['test'], self.workflow.targets)
        test = self.workflow.targets['test']

        self.assertItemsEqual(['foo','bar','baz'], test.input)
        self.assertItemsEqual(['qux','quux'], test.output)