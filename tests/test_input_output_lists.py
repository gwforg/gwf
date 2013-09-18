import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))

def _fname(workflow, fname):
    return os.path.normpath(os.path.join(workflow.working_dir, fname))



class ListTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'input-output-lists.gwf'))

    def test_lists(self):
        self.assertItemsEqual(['infiles','outfiles'], self.workflow.lists)
        
    def test_target(self):
        self.assertItemsEqual(['test','test_new'], self.workflow.targets)
        test = self.workflow.targets['test']

        self.assertItemsEqual([_fname(self.workflow,fn) for fn in ['foo','bar','baz']],
                              test.input)
        self.assertItemsEqual([_fname(self.workflow,fn) for fn in ['qux','quux']],
                              test.output)
                              
        test = self.workflow.targets['test_new']

        self.assertItemsEqual([_fname(self.workflow,fn) for fn in ['foo_new','bar_new','baz_new']],
                              test.input)
        self.assertItemsEqual([_fname(self.workflow,fn) for fn in ['qux_new','quux_new']],
                              test.output)