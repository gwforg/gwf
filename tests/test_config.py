import pytest

from gwf.conf import FileConfig


def test_load_config_from_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write('{"foo": "bar"}\n')

    c = FileConfig.load(config_path)
    assert len(c) == 4
    assert c["foo"] == "bar"


def test_load_config_from_nonexisting_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    c = FileConfig.load(config_path)
    assert len(c) == 3


def test_dump_config_to_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write('{"foo": "bar"}\n')

    c1 = FileConfig.load(config_path)
    c1["baz"] = "foo"
    c1.dump()

    c2 = FileConfig.load(str(config_path))
    assert len(c2) == 5
    assert c2["baz"] == "foo"
    assert c2["foo"] == "bar"
