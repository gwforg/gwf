import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(os.path.realpath(__file__))

def _fname(workflow, fname):
    return os.path.normpath(os.path.join(workflow.working_dir, fname))



class SimpleWorkflowTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'simple_pipe.gwf'))
        
    def test_working_dir(self):
        self.assertTrue(os.path.samefile(self.workflow.working_dir, testdir))
        
    def test_targets(self):
        targets = self.workflow.targets
        assert len(targets) == 3
        
        self.assertIn('SingleSource', targets)
        self.assertIn('SinglePipe', targets)
        self.assertIn('SingleSink', targets)
        
        source = targets['SingleSource']
        self.assertTrue(source.can_execute, "Tasks should be able to execute")
        self.assertTrue(source.should_run, 
                         "Should run if we haven't made the output file")
                         
        self.assertEqual(len(source.input), 0, 
                         "SingleSource has no input files")
        self.assertEqual(len(source.output), 1, 
                         "SingleSource has one output file")
        self.assertEqual(len(source.pbs_options), 0,
                         "SingleSource has no PBS options")
        

        pipe = targets['SinglePipe']
        self.assertTrue(pipe.can_execute, "Tasks should be able to execute")
        self.assertTrue(pipe.should_run, 
                         "Should run if we haven't made the output file")
                         
        self.assertEqual(len(pipe.input), 1,
                         "SinglePipe has one input file")
        self.assertEqual(len(pipe.output), 1, 
                         "SinglePipe has one output file")
        self.assertEqual(len(pipe.pbs_options), 0,
                         "SinglePipe has no PBS options")


        sink = targets['SingleSink']
        self.assertTrue(sink.can_execute, "Tasks should be able to execute")
        self.assertTrue(sink.should_run, 
                         "Should run since there is no output file")
                         
        self.assertEqual(len(sink.input), 1 ,
                         "SingleSink has one input file")
        self.assertEqual(len(sink.output), 0, 
                         "SingleSink has no output files")
        self.assertEqual(len(sink.pbs_options), 0,
                         "SingleSink has no PBS options")

        
    def test_scripting(self):
        # Not sure how to test scripts properly, but at least make sure
        # I can write the scripts in the right place. This should thrown
        # an exception or something if it doesn't work.
        for target in self.workflow.targets.values():
            target.write_script()
            self.assertTrue(target.script_name.startswith(os.path.join(testdir,'.scripts')))

    def test_dependencies(self):
        targets = self.workflow.targets
        
        source = targets['SingleSource']
        self.assertEqual(len(source.dependencies), 0)
        
        pipe = targets['SinglePipe']
        self.assertEqual(len(pipe.dependencies), 1)
        self.assertEqual(pipe.dependencies[0], (_fname(self.workflow,'some_file'), source))
        
        sink = targets['SingleSink']
        self.assertEqual(len(sink.dependencies), 1)
        self.assertEqual(sink.dependencies[0], (_fname(self.workflow,'some_other_file'), pipe))
        
        
        
class SystemFileWorkflowTester(unittest.TestCase):

    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'system_files.gwf'))
        # create one of the two system files
        open(os.path.join(testdir,'sysfile1'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'sysfile1'))
        
    def test_system_files_input(self):
        target = self.workflow.targets['Target']
        self.assertEqual(len(target.dependencies), 2)
        deps_table = dict(target.dependencies)
        self.assertIn(_fname(self.workflow,'sysfile1'), deps_table)
        self.assertIn(_fname(self.workflow,'sysfile2'), deps_table)
        
        sysfile1 = deps_table[_fname(self.workflow,'sysfile1')]
        self.assertTrue(sysfile1.file_exists, 
                        'We created it so it should be there')
        self.assertFalse(sysfile1.should_run,
                         'By convension this is false if the file exists')
        
        sysfile2 = deps_table[_fname(self.workflow,'sysfile2')]
        self.assertFalse(sysfile2.file_exists, 
                         'We did not created it so it should not be there')
        self.assertTrue(sysfile2.should_run,
                        'By convension this is true when the file is missing')
        