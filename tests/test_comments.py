import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class CommentsTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'comments.gwf'))

    def test_lists(self):
        self.assertItemsEqual([], self.workflow.lists)
        
    def test_target(self):
        self.assertItemsEqual([], self.workflow.targets)

