import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))


class TemplateTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'templates.gwf'))
        
    def test_templates(self):
        self.assertEqual(len(self.workflow.templates), 1)
        self.assertIn('cat', self.workflow.templates)
        cat_template = self.workflow.templates['cat']
        self.assertItemsEqual(['file1','file2'], cat_template.parameters)
        
    def test_template_targets(self):
        self.assertEqual(len(self.workflow.template_targets), 1)
        self.assertIn('cat-foo-bar', self.workflow.template_targets)
        tt = self.workflow.template_targets['cat-foo-bar']
        self.assertEqual(tt.template, 'cat')
        self.assertItemsEqual(['file1','file2'], tt.assignments.keys())
        self.assertItemsEqual(['foo','bar'], tt.assignments.values())
        self.assertEqual('foo',tt.assignments['file1'])
        self.assertEqual('bar',tt.assignments['file2'])
        
    def test_targets(self):
        self.assertEqual(len(self.workflow.targets), 1)
        target = self.workflow.targets['cat-foo-bar']
        self.assertItemsEqual(['foo'],target.input)
        self.assertItemsEqual(['bar'],target.output)