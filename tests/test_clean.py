import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))

def _fname(workflow, fname):
    return os.path.normpath(os.path.join(workflow.working_dir, fname))



class ListTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'input-output-lists.gwf'))

    def tearDown(self):
        qux = _fname(self.workflow, 'qux')
        quux = _fname(self.workflow, 'quux')
        if os.path.exists(qux): os.remove(qux)
        if os.path.exists(quux): os.remove(quux)

    def test_target(self):
        self.assertItemsEqual(['test','test_new'], self.workflow.targets)
        test = self.workflow.targets['test']
        self.assertItemsEqual([_fname(self.workflow,fn) 
                                for fn in ['qux','quux']],
                              test.output)

        test = self.workflow.targets['test_new']
        self.assertItemsEqual([_fname(self.workflow,fn) 
                                for fn in ['qux_new','quux_new']],
                              test.output)
                              
    def test_existing_output(self):
        self.assertItemsEqual(['test','test_new'], self.workflow.targets)
        test = self.workflow.targets['test']

        existing = test.get_existing_outfiles()
        expected = []
        self.assertItemsEqual(existing, expected)

        qux = _fname(self.workflow, 'qux')
        quux = _fname(self.workflow, 'quux')
        
        open(qux,'w').close()
        existing = test.get_existing_outfiles()
        expected = [qux]
        self.assertItemsEqual(existing, expected)

        open(quux,'w').close()
        existing = test.get_existing_outfiles()
        expected = [qux,quux]
        self.assertItemsEqual(existing, expected)

        test.clean_target()
        existing = test.get_existing_outfiles()
        expected = []
        self.assertItemsEqual(existing, expected)
        
        self.tearDown()

    def test_existing_output_new(self):
        self.assertItemsEqual(['test','test_new'], self.workflow.targets)
        test = self.workflow.targets['test_new']

        existing = test.get_existing_outfiles()
        expected = []
        self.assertItemsEqual(existing, expected)

        qux = _fname(self.workflow, 'qux_new')
        quux = _fname(self.workflow, 'quux_new')
        
        open(qux,'w').close()
        existing = test.get_existing_outfiles()
        expected = [qux]
        self.assertItemsEqual(existing, expected)

        open(quux,'w').close()
        existing = test.get_existing_outfiles()
        expected = [qux,quux]
        self.assertItemsEqual(existing, expected)

        test.clean_target()
        existing = test.get_existing_outfiles()
        expected = []
        self.assertItemsEqual(existing, expected)

        self.tearDown()