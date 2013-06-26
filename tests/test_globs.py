import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class GlobTester(unittest.TestCase):

    def setUp(self):
        open(os.path.join(testdir,'glob.1'),'a').close()
        open(os.path.join(testdir,'glob.2'),'a').close()
        open(os.path.join(testdir,'glob.3'),'a').close()
        self.workflow = parse(os.path.join(testdir,'globs.gwf'))

    def tearDown(self):
        os.remove(os.path.join(testdir,'glob.1'))
        os.remove(os.path.join(testdir,'glob.2'))
        os.remove(os.path.join(testdir,'glob.3'))
        
    def test_lists(self):
        self.assertItemsEqual(['globs'], self.workflow.lists)
        globs = self.workflow.lists['globs']
        files = [os.path.join(testdir,f) for f in  ['glob.1','glob.2','glob.3']]
        self.assertItemsEqual(files, globs.elements)
        
