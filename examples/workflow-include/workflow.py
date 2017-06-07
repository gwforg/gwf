from gwf import Workflow

gwf = Workflow()
gwf.include('other_workflow/workflow.py', namespace='other')
gwf.target('World', inputs=['other_workflow/a.txt'], outputs=['b.txt']) << """
cat other_workflow/a.txt > b.txt
echo world >> b.txt
"""

