from gwf.conf import FileConfig


def test_load_config_from_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write('{"foo": "bar"}\n')

    c = FileConfig.load(config_path)
    assert len(c) == 5
    assert c["foo"] == "bar"
    assert c["require_all_outputs"] is True


def test_load_config_from_nonexisting_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    c = FileConfig.load(config_path)
    assert len(c) == 4
    assert c["require_all_outputs"] is True


def test_dump_config_to_file(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write('{"foo": "bar"}\n')

    c1 = FileConfig.load(config_path)
    c1["baz"] = "foo"
    c1.dump()

    c2 = FileConfig.load(str(config_path))
    assert len(c2) == 6
    assert c2["baz"] == "foo"
    assert c2["foo"] == "bar"


def test_unset_overridden_default_restores_default(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write('{"clean_logs": false}\n')

    c = FileConfig.load(config_path)
    del c["clean_logs"]

    assert c["clean_logs"] is True


def test_unset_key_with_only_a_default_keeps_default(tmpdir):
    c = FileConfig.load(tmpdir.join(".gwfconf.json"))

    del c["clean_logs"]

    assert c["clean_logs"] is True


def test_get_namespace_only_matches_whole_namespace(tmpdir):
    config_path = tmpdir.join(".gwfconf.json")
    config_path.write(
        '{"backend.slurm.log_mode": "merged", '
        '"backend.slurmish.log_mode": "none", '
        '"backend.slurm": "oops"}\n'
    )

    c = FileConfig.load(config_path)

    assert c.get_namespace("backend.slurm") == {"log_mode": "merged"}
