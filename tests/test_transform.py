import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class ShellTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'transform.gwf'))

        
    def test_lists(self):
        self.assertItemsEqual(['foos','foos2','bars','bars2'],
                              self.workflow.lists)
                              
        foos = self.workflow.lists['foos']
        bars = self.workflow.lists['bars']
        self.assertItemsEqual(['foo1','foo2','foo3'], foos.elements)
        self.assertItemsEqual(['bar1','bar2','bar3'], bars.elements)
        
        foos = self.workflow.lists['foos2']
        bars = self.workflow.lists['bars2']
        self.assertItemsEqual(['foo1','foo2','foo3'], foos.elements)
        self.assertItemsEqual(['bar1','bar2','bar3'], bars.elements)
