import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))

def _fname(workflow, fname):
    return os.path.normpath(os.path.join(workflow.working_dir, fname))



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
        self.assertIn('cat_foo_bar', self.workflow.template_targets)
        tt = self.workflow.template_targets['cat_foo_bar']
        self.assertEqual(tt.template, 'cat')
        self.assertItemsEqual(['file1','file2'], tt.assignments.keys())
        self.assertItemsEqual(['foo','bar'], tt.assignments.values())
        self.assertEqual('foo',tt.assignments['file1'])
        self.assertEqual('bar',tt.assignments['file2'])
        
    def test_targets(self):
        targets = self.workflow.targets
        self.assertEqual(len(targets), 3)
        self.assertItemsEqual(['cat_foo_bar','make_foo','eat_bar'], targets)
        
        template_target = targets['cat_foo_bar']
        self.assertItemsEqual([_fname(self.workflow,'foo')],template_target.input)
        self.assertItemsEqual([_fname(self.workflow,'bar')],template_target.output)
        
        foo = targets['make_foo']
        bar = targets['eat_bar']
        
        self.assertEqual(len(template_target.dependencies), 1)
        self.assertEqual(template_target.dependencies[0], (_fname(self.workflow,'foo'),foo))
        
        self.assertEqual(len(bar.dependencies), 1)
        self.assertEqual(bar.dependencies[0], (_fname(self.workflow,'bar'),template_target))
        