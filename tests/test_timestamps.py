import unittest
import os, os.path

from gwf.parser import parse
testdir = os.path.dirname(__file__)

class SourceMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
    
    def test_source_should_run(self):
        source = self.workflow.targets['source']
        self.assertTrue(source.should_run,
                        "Output file is missing, so the target should run.")
        

class SourceExisting(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'source'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'source'))
    
    def test_source_should_run(self):
        source = self.workflow.targets['source']
        self.assertFalse(source.should_run,
                         "Output file exists so the target should not run.")
                         

                         
class SinkMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
    def test_sink_should_run(self):
        sink = self.workflow.targets['sink']
        self.assertTrue(sink.should_run, "A sink should always run!")
        

class SinkExisting(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'sink'),'a').close()
    def tearDown(self):
        os.remove(os.path.join(testdir,'sink'))
    def test_sink_should_run(self):
        sink = self.workflow.targets['sink']
        self.assertTrue(sink.should_run, "A sink should always run!")


class PipeBothMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
    def tearDown(self):
        pass
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeInMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'file1'),'a').close()
    def tearDown(self):
        os.remove(os.path.join(testdir,'file1'))
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeOutMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'file2'),'a').close()
    def tearDown(self):
        os.remove(os.path.join(testdir,'file2'))
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeInYoungest(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'file2'),'a').close()
        os.system('sleep 1') # make sure they have different time stamps
        open(os.path.join(testdir,'file1'),'a').close()
    def tearDown(self):
        os.remove(os.path.join(testdir,'file1'))
        os.remove(os.path.join(testdir,'file2'))
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, 
                        "The infile is younger than the outfile "
                        "so we should run!")

class PipeOutYoungest(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'file1'),'a').close()
        os.system('sleep 1') # make sure they have different time stamps
        open(os.path.join(testdir,'file2'),'a').close()
    def tearDown(self):
        os.remove(os.path.join(testdir,'file1'))
        os.remove(os.path.join(testdir,'file2'))
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertFalse(pipe.should_run, 
                         "The outfile is older than the infile "
                         "so we should not run!")



class TwoIn_12o(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'in1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in2'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in1'))
        os.remove(os.path.join(testdir,'in2'))
        os.remove(os.path.join(testdir,'out'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertFalse(target.should_run, target.reason_to_run)

class TwoIn_1o2(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'in1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in2'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in1'))
        os.remove(os.path.join(testdir,'in2'))
        os.remove(os.path.join(testdir,'out'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertTrue(target.should_run, target.reason_to_run)

class TwoIn_o12(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'out'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in2'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in1'))
        os.remove(os.path.join(testdir,'in2'))
        os.remove(os.path.join(testdir,'out'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertTrue(target.should_run, target.reason_to_run)



class TwoOut_i12(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'in'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out2'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in'))
        os.remove(os.path.join(testdir,'out1'))
        os.remove(os.path.join(testdir,'out2'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertFalse(target.should_run, target.reason_to_run)

class TwoOut_1i2(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'out1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out2'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in'))
        os.remove(os.path.join(testdir,'out1'))
        os.remove(os.path.join(testdir,'out2'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertTrue(target.should_run, target.reason_to_run)

class TwoOut_12i(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(os.path.join(testdir,'timestamps.gwf'))
        open(os.path.join(testdir,'out1'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'out2'),'a').close()
        os.system('sleep 1')
        open(os.path.join(testdir,'in'),'a').close()
        
    def tearDown(self):
        os.remove(os.path.join(testdir,'in'))
        os.remove(os.path.join(testdir,'out1'))
        os.remove(os.path.join(testdir,'out2'))
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertTrue(target.should_run, target.reason_to_run)
