from gwf import Workflow

gwf = Workflow()

gwf.target('SayHello', inputs=['name.txt'], outputs=['greeting.txt']) << """
echo -n "Hello " > greeting.txt
cat name.txt >> greeting.txt
"""
