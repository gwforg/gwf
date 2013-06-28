import unittest
import os, os.path

from gwf.parser import parse
testdir = os.path.dirname(__file__)

def _fname(workflow, fname):
    return os.path.normpath(os.path.join(workflow.working_dir, fname))


class DependencyGraphTest(unittest.TestCase):

    def setUp(self):
        open(os.path.join(testdir,'e21'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e22'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e31'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e32'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'final'),'a').close()
        self.workflow = parse(os.path.join(testdir,'graph.gwf'))
        self.graph = self.workflow.dependency_graph
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'e21'))
        os.remove(os.path.join(testdir,'e22'))
        os.remove(os.path.join(testdir,'e31'))
        os.remove(os.path.join(testdir,'e32'))
        os.remove(os.path.join(testdir,'final'))
        
        
    def test_nodes(self):
        nodes = self.graph.nodes
        targets = self.workflow.targets
        
        self.assertEqual(len(nodes), 5)
        
        n1 = nodes['n1']
        self.assertEqual(len(n1.dependencies), 0)
        self.assertEqual(n1.name, 'n1')
        self.assertEqual(n1.task, targets['n1'])
        self.assertTrue(n1.task.should_run)
        self.assertTrue(n1.should_run)
        
        n2 = nodes['n2']
        self.assertEqual(len(n2.dependencies), 1)
        self.assertEqual(n2.dependencies[0], (_fname(self.workflow,'e1'),n1))
        self.assertEqual(n2.name, 'n2')
        self.assertEqual(n2.task, targets['n2'])
        self.assertTrue(n2.task.should_run)
        self.assertTrue(n2.should_run)
        
        n31 = nodes['n31']
        self.assertEqual(len(n31.dependencies), 1)
        self.assertEqual(n31.dependencies[0], (_fname(self.workflow,'e21'),n2))
        self.assertEqual(n31.name, 'n31')
        self.assertEqual(n31.task, targets['n31'])
        self.assertFalse(n31.task.should_run, n31.task.reason_to_run)
        self.assertTrue(n31.should_run)
        
        n32 = nodes['n32']
        self.assertEqual(len(n32.dependencies), 1)
        self.assertEqual(n32.dependencies[0], (_fname(self.workflow,'e22'),n2))
        self.assertEqual(n32.name, 'n32')
        self.assertEqual(n32.task, targets['n32'])
        self.assertFalse(n32.task.should_run, n32.task.reason_to_run)
        self.assertTrue(n32.should_run)
        
        n4 = nodes['n4']
        self.assertEqual(len(n4.dependencies), 2)
        deps = dict(n4.dependencies)
        self.assertEqual(deps[_fname(self.workflow,'e31')], n31)
        self.assertEqual(deps[_fname(self.workflow,'e32')], n32)
        self.assertEqual(n4.name, 'n4')
        self.assertEqual(n4.task, targets['n4'])
        self.assertFalse(n4.task.should_run, n4.task.reason_to_run)
        self.assertTrue(n4.should_run)
        
class ScheduleTest1(unittest.TestCase):

    def setUp(self):
        open(os.path.join(testdir,'e21'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e22'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e31'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e32'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'final'),'a').close()
        self.workflow = parse(os.path.join(testdir,'graph.gwf'))
        self.graph = self.workflow.dependency_graph
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'e21'))
        os.remove(os.path.join(testdir,'e22'))
        os.remove(os.path.join(testdir,'e31'))
        os.remove(os.path.join(testdir,'e32'))
        os.remove(os.path.join(testdir,'final'))
        
    def test_schedule(self):
        schedule, scheduled = self.graph.schedule('n4')
        self.assertEqual(len(schedule), 5)
        self.assertItemsEqual(['n1','n2','n31','n32','n4'], scheduled)
        
        n1 = self.graph.nodes['n1']
        n2 = self.graph.nodes['n2']
        n31 = self.graph.nodes['n31']
        n32 = self.graph.nodes['n32']
        n4 = self.graph.nodes['n4']
        self.assertItemsEqual([n1,n2,n31,n32,n4], schedule)
        self.assertLess(schedule.index(n1),schedule.index(n2))
        self.assertLess(schedule.index(n2),schedule.index(n31))
        self.assertLess(schedule.index(n2),schedule.index(n32))
        self.assertLess(schedule.index(n31),schedule.index(n4))
        self.assertLess(schedule.index(n32),schedule.index(n4))

class ScheduleTest2(unittest.TestCase):

    def setUp(self):
        open(os.path.join(testdir,'e1'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e22'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e31'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e32'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'final'),'a').close()
        self.workflow = parse(os.path.join(testdir,'graph.gwf'))
        self.graph = self.workflow.dependency_graph
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'e1'))
        os.remove(os.path.join(testdir,'e22'))
        os.remove(os.path.join(testdir,'e31'))
        os.remove(os.path.join(testdir,'e32'))
        os.remove(os.path.join(testdir,'final'))
        
    def test_schedule(self):
        schedule, scheduled = self.graph.schedule('n4')
        self.assertEqual(len(schedule), 4)
        self.assertItemsEqual(['n2','n31','n32','n4'], scheduled)
        
        n2 = self.graph.nodes['n2']
        n31 = self.graph.nodes['n31']
        n32 = self.graph.nodes['n32']
        n4 = self.graph.nodes['n4']
        self.assertItemsEqual([n2,n31,n32,n4], schedule)
        self.assertLess(schedule.index(n2),schedule.index(n31))
        self.assertLess(schedule.index(n2),schedule.index(n32))
        self.assertLess(schedule.index(n31),schedule.index(n4))
        self.assertLess(schedule.index(n32),schedule.index(n4))


class ScheduleTest3(unittest.TestCase):

    def setUp(self):
        open(os.path.join(testdir,'e1'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e21'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e22'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e31'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'e32'),'a').close() ; os.system('sleep 1')
        open(os.path.join(testdir,'final'),'a').close()
        self.workflow = parse(os.path.join(testdir,'graph.gwf'))
        self.graph = self.workflow.dependency_graph
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'e1'))
        os.remove(os.path.join(testdir,'e21'))
        os.remove(os.path.join(testdir,'e22'))
        os.remove(os.path.join(testdir,'e31'))
        os.remove(os.path.join(testdir,'e32'))
        os.remove(os.path.join(testdir,'final'))
        
    def test_schedule(self):
        schedule, scheduled = self.graph.schedule('n4')
        self.assertEqual(len(schedule), 0)
