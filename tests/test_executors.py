import io
from pathlib import Path

import pytest

from gwf.core import Target
from gwf.exceptions import GWFError
from gwf.executors import Conda, Pixi, deserialize, serialize


def test_pixi_project_path_must_be_absolute(tmp_path):
    project = tmp_path / "project"

    assert Pixi(project=project).project == project

    with pytest.raises(GWFError, match="Pixi project path must be absolute"):
        Pixi(project=Path("project"))


def test_serialize_describes_target_in_comment_block():
    target = Target(
        "Align",
        inputs={
            "reads": ["raw/sample1_R1.fastq.gz", "raw/sample1_R2.fastq.gz"],
            "ref": "ref/genome.fa",
        },
        outputs=["aligned.bam"],
        options={"cores": 4},
        working_dir="/data/project",
        executor=Conda(env="bwa"),
        spec="bwa mem ref/genome.fa raw/*.fastq.gz > aligned.bam",
    )
    buffer = io.StringIO()
    serialize(target, buffer)

    comment = [
        line.removeprefix("#GWF COMMENT ")
        for line in buffer.getvalue().splitlines()
        if line.startswith("#GWF COMMENT ")
    ]
    assert comment == [
        "Target(",
        "    name='Align',",
        "    inputs={'reads': ['raw/sample1_R1.fastq.gz', 'raw/sample1_R2.fastq.gz'],",
        "            'ref': 'ref/genome.fa'},",
        "    outputs=['aligned.bam'],",
        "    options={'cores': 4},",
        "    working_dir='/data/project',",
        "    executor=Conda(env='bwa', debug_mode=False),",
        ")",
    ]

    buffer.seek(0)
    restored = deserialize(buffer)
    assert restored.name == target.name
    assert restored.spec == target.spec
