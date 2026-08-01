from pathlib import Path

import pytest

from gwf.exceptions import GWFError
from gwf.executors import Pixi


def test_pixi_project_path_must_be_absolute(tmp_path):
    project = tmp_path / "project"

    assert Pixi(project=project).project == project

    with pytest.raises(GWFError, match="Pixi project path must be absolute"):
        Pixi(project=Path("project"))
