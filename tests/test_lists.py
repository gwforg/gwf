import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class ListTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'lists.gwf'))

    def test_lists(self):
        self.assertItemsEqual(['singleton','list1','list2','list3'],
                              self.workflow.lists)
                              
        singleton = self.workflow.lists['singleton']
        self.assertEqual('singleton', singleton.name)
        self.assertItemsEqual(['x'], singleton.elements)
        
        list1 = self.workflow.lists['list1']
        self.assertEqual('list1', list1.name)
        self.assertItemsEqual(['elm1','elm2','elm3'], list1.elements)
        
        list2 = self.workflow.lists['list2']
        self.assertEqual('list2', list2.name)
        self.assertItemsEqual(['x','y','z'], list2.elements)
        
    def test_targets(self):
        self.assertItemsEqual(['singleton',
                                'one-elm1','one-elm2','one-elm3',
                                'two-single-a','two-single-b','two-single-c',
                                'two-elm1-x','two-elm2-y','two-elm3-z'],
                              self.workflow.targets)
                              
        one = self.workflow.targets['one-elm1']
        self.assertItemsEqual(['elm1'],one.input)
        self.assertItemsEqual(['elm1.out'],one.output)

        one = self.workflow.targets['one-elm2']
        self.assertItemsEqual(['elm2'],one.input)
        self.assertItemsEqual(['elm2.out'],one.output)

        one = self.workflow.targets['one-elm3']
        self.assertItemsEqual(['elm3'],one.input)
        self.assertItemsEqual(['elm3.out'],one.output)

        
        two = self.workflow.targets['two-elm1-x']
        self.assertItemsEqual(['elm1'],two.input)
        self.assertItemsEqual(['x'],two.output)
        
        two = self.workflow.targets['two-elm2-y']
        self.assertItemsEqual(['elm2'],two.input)
        self.assertItemsEqual(['y'],two.output)
        
        two = self.workflow.targets['two-elm3-z']
        self.assertItemsEqual(['elm3'],two.input)
        self.assertItemsEqual(['z'],two.output)
        
