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
        targets = self.workflow.targets
        self.assertEqual(len(targets), 3)
        self.assertItemsEqual(['cat-foo-bar','make-foo','eat-bar'], targets)
        
        template_target = targets['cat-foo-bar']
        self.assertItemsEqual(['foo'],template_target.input)
        self.assertItemsEqual(['bar'],template_target.output)
        
        foo = targets['make-foo']
        bar = targets['eat-bar']
        
        self.assertEqual(len(template_target.dependencies), 1)
        self.assertEqual(template_target.dependencies[0], ('foo',foo))
        
        self.assertEqual(len(bar.dependencies), 1)
        self.assertEqual(bar.dependencies[0], ('bar',template_target))
        