import pytest

from gwf.utils import parse_path, cache, PersistableDict, ensure_trailing_newline, retry


def test_cache_returns_same_object_when_called_twice_with_same_args():
    func = cache(lambda x: object())
    obj1 = func(42)
    obj2 = func(42)
    assert id(obj1) == id(obj2)


def test_cache_does_not_return_same_object_when_called_with_diff_args():
    func = cache(lambda x: object())
    obj1 = func(42)
    obj2 = func(43)
    assert id(obj1) != id(obj2)


@pytest.mark.parametrize(
    "path,parsed_path",
    [
        ("/some/dir/workflow.py", ("/some/dir", "workflow", "gwf")),
        ("/some/dir/workflow.py:other", ("/some/dir", "workflow", "other")),
        ("/some/dir/other.py:other", ("/some/dir", "other", "other")),
        (":other", ("", "workflow", "other")),
    ],
)
def test_parse_path(path, parsed_path):
    assert parse_path(path) == parsed_path


@pytest.mark.parametrize(
    "default_obj,default_file,parsed_path",
    [
        ("wf1", "file1.py", ("", "file1", "wf1")),
        ("wf1", "file2.py", ("", "file2", "wf1")),
        ("wf2", "file1.py", ("", "file1", "wf2")),
        ("wf2", "file2.py", ("", "file2", "wf2")),
    ],
)
def test_parse_path_defaults(default_obj, default_file, parsed_path):
    assert (
        parse_path(":", default_obj=default_obj, default_file=default_file)
        == parsed_path
    )


def test_dict_is_empty_if_file_does_not_exist(tmpdir):
    with tmpdir.as_cwd():
        d = PersistableDict("test.json")
        assert d == {}


def test_dict_write_read(tmpdir):
    with tmpdir.as_cwd():
        d1 = PersistableDict("test.json")
        d1["foo"] = "bar"
        d1.persist()

        assert tmpdir.join("test.json").exists()

        d2 = PersistableDict("test.json")
        assert d2 == {"foo": "bar"}


def test_ensure_trailing_newline():
    assert ensure_trailing_newline("") == "\n"
    assert ensure_trailing_newline("foo\nbar\n") == "foo\nbar\n"
    assert ensure_trailing_newline("foo\nbar") == "foo\nbar\n"


def test_retry_1(mocker, no_sleep):
    failing_func = mocker.Mock(side_effect=[Exception, Exception, 42])
    failing_func.__name__ = "failing_func"

    retry_func = retry(on_exc=Exception, max_retries=5)
    wrapped_func = retry_func(failing_func)

    assert wrapped_func() == 42
    assert len(failing_func.call_args_list) == 3


def test_retry_2(mocker, no_sleep):
    failing_func = mocker.Mock(side_effect=[Exception, Exception, Exception])
    failing_func.__name__ = "failing_func"

    retry_func = retry(on_exc=Exception, max_retries=3)
    wrapped_func = retry_func(failing_func)

    with pytest.raises(retry.RetryError):
        wrapped_func()

    assert len(failing_func.call_args_list) == 3


def test_retry_3(mocker, no_sleep):
    failing_func = mocker.Mock(side_effect=[Exception, Exception, 42])
    failing_func.__name__ = "failing_func"

    callback_func = mocker.Mock()

    retry_func = retry(on_exc=Exception, max_retries=3, callback=callback_func)
    wrapped_func = retry_func(failing_func)

    assert wrapped_func() == 42
    assert len(failing_func.call_args_list) == 3
    assert len(callback_func.call_args_list) == 2


def test_retry_4(mocker, no_sleep):
    succeeding_func = mocker.Mock(side_effect=[42])
    succeeding_func.__name__ = "succeeding_func"

    retry_func = retry(on_exc=Exception, max_retries=5)
    wrapped_func = retry_func(succeeding_func)

    assert wrapped_func() == 42
    assert len(succeeding_func.call_args_list) == 1
