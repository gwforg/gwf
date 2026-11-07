import pytest

from gwf.utils import ensure_trailing_newline, is_valid_name


def test_ensure_trailing_newline():
    assert ensure_trailing_newline("") == "\n"
    assert ensure_trailing_newline("foo\nbar\n") == "foo\nbar\n"
    assert ensure_trailing_newline("foo\nbar") == "foo\nbar\n"


@pytest.mark.parametrize("name", ["foo", "_foo", "Foo_1", "sample.1"])
def test_is_valid_name_accepts_valid_names(name):
    assert is_valid_name(name)


@pytest.mark.parametrize("name", ["", "1foo", "foo bar", "foo-bar", "foo\n", "\nfoo"])
def test_is_valid_name_rejects_invalid_names(name):
    assert not is_valid_name(name)
