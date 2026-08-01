from gwf import IsolationConfig, Workflow


gwf = Workflow()


copy_spec = """
test ! -L inputs/samples/sample.txt
cat inputs/samples/sample.txt > results/copied/output.txt
mkdir -p work/copied
echo "This file only exists inside the isolated directory." > work/copied/undeclared.txt
"""
copy_target = gwf.target(
    "CopyInput",
    inputs=["inputs/samples/sample.txt"],
    outputs=["results/copied/output.txt"],
    isolation=IsolationConfig(inputs="copy"),
)
copy_target << copy_spec


symlink_spec = """
test -L inputs/references/reference.txt
cat inputs/references/reference.txt > results/symlinked/output.txt
mkdir -p work/symlinked
echo "This file only exists inside the isolated directory." > work/symlinked/undeclared.txt
"""
symlink_target = gwf.target(
    "SymlinkInput",
    inputs=["inputs/references/reference.txt"],
    outputs=["results/symlinked/output.txt"],
    isolation=IsolationConfig(inputs="symlink"),
)
symlink_target << symlink_spec


nested_spec = """
test ! -L inputs/samples/sample.txt
test ! -L inputs/references/reference.txt
cat inputs/samples/sample.txt inputs/references/reference.txt > results/nested/combined.txt
mkdir -p work/nested
cat inputs/samples/sample.txt > work/nested/intermediate.txt
"""
nested_target = gwf.target(
    "NestedHierarchy",
    inputs={
        "sample": "inputs/samples/sample.txt",
        "reference": "inputs/references/reference.txt",
    },
    outputs={"combined": "results/nested/combined.txt"},
    isolation=IsolationConfig(inputs="copy"),
)
nested_target << nested_spec


missing_output_spec = """
mkdir -p work/missing-output
cat inputs/samples/sample.txt > work/missing-output/produced-but-undeclared.txt
"""
missing_output_target = gwf.target(
    "MissingOutput",
    inputs=["inputs/samples/sample.txt"],
    outputs=["results/missing-output.txt"],
    isolation=True,
)
missing_output_target << missing_output_spec


missing_input_spec = """
echo "This command is never run." > missing-input-command-ran.txt
"""
missing_input_target = gwf.target(
    "MissingInput",
    inputs=["inputs/samples/does-not-exist.txt"],
    outputs=["results/missing-input.txt"],
    isolation=True,
)
missing_input_target << missing_input_spec

# Only declared outputs from successful targets are moved back under results/.
# Files created under work/ are discarded with the temporary directories.
# MissingOutput fails because it does not create results/missing-output.txt.
# MissingInput is skipped before execution because its declared input is absent.
