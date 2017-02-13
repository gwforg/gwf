from gwf import Workflow
from explicit import *

gwf = Workflow()

gwf // target('SayHello') // inputs('name.txt') // outputs('greeting.txt') \
    << """
echo -n "Hello " > greeting.txt
cat name.txt >> greeting.txt
"""

gwf // target("World") // inputs('greeting.txt') // outputs('world.txt') \
    << """
cat greeting.txt > world.txt
echo "world" >> world.txt
"""

gwf // target("Universe") // inputs('greeting.txt') // outputs('universe.txt') \
    << """
cat greeting.txt > universe.txt
echo "universe" >> universe.txt
"""

gwf // target("All") \
    // inputs('world.txt', 'universe.txt') \
    // outputs('all.txt') \
    << """
cat world.txt > all.txt
cat universe.txt >> all.txt
"""
