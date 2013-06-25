
import os
import sys
from gwf.parser import parse

workflow = parse(sys.argv[1])
print workflow.templates