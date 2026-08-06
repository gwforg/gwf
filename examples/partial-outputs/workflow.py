"""Demonstrate scheduling a dependency from only the outputs it needs."""

from gwf import Workflow


gwf = Workflow()


gwf.target(
    "ProduceInputs",
    inputs=[],
    outputs=["x.txt", "y.txt"],
) << """
printf 'data for x\\n' > x.txt
printf 'data for y\\n' > y.txt
"""

gwf.target(
    "UseX",
    inputs=["x.txt"],
    outputs=["x-result.txt"],
) << """
cat x.txt > x-result.txt
"""

gwf.target(
    "UseY",
    inputs=["y.txt"],
    outputs=["y-result.txt"],
) << """
cat y.txt > y-result.txt
"""
