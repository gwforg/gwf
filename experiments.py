import os
import sys
from gwf.parser import parse

workflow = parse(sys.argv[1])

print workflow.lists
print workflow.templates
print workflow.template_targets
print workflow.targets

