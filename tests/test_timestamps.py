import unittest
import os 

from gwf.parser import parse
testdir = os.path.dirname(__file__)

class SourceMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
    
    def test_source_should_run(self):
        source = self.workflow.targets['source']
        self.assertTrue(source.should_run,
                        "Output file is missing, so the target should run.")
        

class SourceExisting(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/source','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/source')
    
    def test_source_should_run(self):
        source = self.workflow.targets['source']
        self.assertFalse(source.should_run,
                         "Output file exists so the target should not run.")
                         

                         
class SinkMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
    def test_sink_should_run(self):
        sink = self.workflow.targets['sink']
        self.assertTrue(sink.should_run, "A sink should always run!")
        

class SinkExisting(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/sink','a').close()
    def tearDown(self):
        os.remove(testdir+'/sink')
    def test_sink_should_run(self):
        sink = self.workflow.targets['sink']
        self.assertTrue(sink.should_run, "A sink should always run!")


class PipeBothMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
    def tearDown(self):
        pass
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeInMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/file1','a').close()
    def tearDown(self):
        os.remove(testdir+'/file1')
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeOutMissing(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/file2','a').close()
    def tearDown(self):
        os.remove(testdir+'/file2')
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, "Files are missing so we should run!")

class PipeInYoungest(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/file2','a').close()
        os.system('sleep 1') # make sure they have different time stamps
        open(testdir+'/file1','a').close()
    def tearDown(self):
        os.remove(testdir+'/file1')
        os.remove(testdir+'/file2')
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertTrue(pipe.should_run, 
                        "The infile is younger than the outfile "
                        "so we should run!")

class PipeOutYoungest(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/file1','a').close()
        os.system('sleep 1') # make sure they have different time stamps
        open(testdir+'/file2','a').close()
    def tearDown(self):
        os.remove(testdir+'/file1')
        os.remove(testdir+'/file2')
    def test_pipe_should_run(self):
        pipe = self.workflow.targets['pipe']
        self.assertFalse(pipe.should_run, 
                         "The outfile is older than the infile "
                         "so we should not run!")



class TwoIn_12o(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/in1','a').close()
        os.system('sleep 1')
        open(testdir+'/in2','a').close()
        os.system('sleep 1')
        open(testdir+'/out','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in1')
        os.remove(testdir+'/in2')
        os.remove(testdir+'/out')
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertFalse(target.should_run, target.reason_to_run)

class TwoIn_1o2(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/in1','a').close()
        os.system('sleep 1')
        open(testdir+'/out','a').close()
        os.system('sleep 1')
        open(testdir+'/in2','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in1')
        os.remove(testdir+'/in2')
        os.remove(testdir+'/out')
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertTrue(target.should_run, target.reason_to_run)

class TwoIn_o12(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/out','a').close()
        os.system('sleep 1')
        open(testdir+'/in1','a').close()
        os.system('sleep 1')
        open(testdir+'/in2','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in1')
        os.remove(testdir+'/in2')
        os.remove(testdir+'/out')
        
    def test_should_run(self):
        target = self.workflow.targets['twoin']
        self.assertTrue(target.should_run, target.reason_to_run)



class TwoOut_i12(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/in','a').close()
        os.system('sleep 1')
        open(testdir+'/out1','a').close()
        os.system('sleep 1')
        open(testdir+'/out2','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in')
        os.remove(testdir+'/out1')
        os.remove(testdir+'/out2')
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertFalse(target.should_run, target.reason_to_run)

class TwoOut_1i2(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/out1','a').close()
        os.system('sleep 1')
        open(testdir+'/in','a').close()
        os.system('sleep 1')
        open(testdir+'/out2','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in')
        os.remove(testdir+'/out1')
        os.remove(testdir+'/out2')
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertTrue(target.should_run, target.reason_to_run)

class TwoOut_12i(unittest.TestCase):
    def setUp(self):
        self.workflow = parse(testdir+'/timestamps.gwf')
        open(testdir+'/out1','a').close()
        os.system('sleep 1')
        open(testdir+'/out2','a').close()
        os.system('sleep 1')
        open(testdir+'/in','a').close()
        
    def tearDown(self):
        os.remove(testdir+'/in')
        os.remove(testdir+'/out1')
        os.remove(testdir+'/out2')
        
    def test_should_run(self):
        target = self.workflow.targets['twoout']
        self.assertTrue(target.should_run, target.reason_to_run)
