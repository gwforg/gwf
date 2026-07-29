import os
import subprocess
import sys

import pytest

from gwf.backends import lsf, pbs, sge, slurm
from gwf.core import Target


@pytest.mark.parametrize(
    ("ops_factory", "target_defaults"),
    [
        (
            lambda tmp_path: slurm.SlurmOps(
                str(tmp_path), "full", True, slurm.TARGET_DEFAULTS
            ),
            slurm.TARGET_DEFAULTS,
        ),
        (
            lambda tmp_path: sge.SGEOps(str(tmp_path), sge.TARGET_DEFAULTS),
            sge.TARGET_DEFAULTS,
        ),
        (
            lambda tmp_path: lsf.LSFOps(str(tmp_path), lsf.TARGET_DEFAULTS),
            lsf.TARGET_DEFAULTS,
        ),
        (
            lambda tmp_path: pbs.PBSOps(str(tmp_path), pbs.TARGET_DEFAULTS),
            pbs.TARGET_DEFAULTS,
        ),
    ],
)
def test_backend_scripts_execute_targets_with_gwf_exec(
    tmp_path, ops_factory, target_defaults
):
    target = Target(
        "Target",
        inputs=[],
        outputs=["output.txt"],
        options=dict(target_defaults),
        working_dir=str(tmp_path),
        spec="echo completed > output.txt",
    )
    script = ops_factory(tmp_path).compile_script(target)
    script_path = tmp_path / "job-script"
    script_path.write_text(script)
    script_path.chmod(0o700)
    environ = os.environ.copy()
    environ["GWF_EXEC_WORKFLOW_ROOT"] = str(tmp_path)

    result = subprocess.run(
        [script_path],
        cwd=tmp_path,
        env=environ,
        text=True,
        capture_output=True,
        check=False,
    )

    assert script.startswith(f"#!{sys.executable} -mgwf.exec\n")
    assert "#GWF SPEC\n\necho completed > output.txt\n" in script
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "output.txt").read_text() == "completed\n"
