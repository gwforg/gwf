from gwf.utils import ensure_trailing_newline


def test_ensure_trailing_newline():
    assert ensure_trailing_newline("") == "\n"
    assert ensure_trailing_newline("foo\nbar\n") == "foo\nbar\n"
    assert ensure_trailing_newline("foo\nbar") == "foo\nbar\n"
