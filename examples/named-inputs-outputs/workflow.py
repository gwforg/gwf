from gwf import Workflow


gwf = Workflow()

producer = gwf.target(
    name="Producer",
    inputs=[],
    outputs={"A": ["a1", "a2"], "B": "b"},
) << """
    touch a1
    touch a2
    touch b
"""

gwf.target(
    name="Consumer1",
    inputs=producer.outputs["A"],
    outputs=["as"],
) << """
    cat a1 a2 > as
"""

gwf.target(
    name="Consumer2",
    inputs=producer.outputs["B"],
    outputs=["c"],
) << """
    touch c
"""
