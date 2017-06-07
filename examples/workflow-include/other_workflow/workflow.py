from gwf import Workflow

gwf = Workflow()
gwf.target('World', inputs=[], outputs=['a.txt']) << """
echo hello > a.txt
"""
