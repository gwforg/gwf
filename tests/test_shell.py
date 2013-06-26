import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class ShellTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'shell.gwf'))

        
    def test_lists(self):
        self.assertItemsEqual(['list1'], self.workflow.lists)
        list1 = self.workflow.lists['list1']
        self.assertItemsEqual(map(str,range(1,6)), list1.elements)
        
