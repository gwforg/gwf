import pytest

from gwf.cli import main
from gwf.plugins.config import humanbool, cast_value


def test_set_get_unset(cli_runner):
    cli_runner.invoke(main, ['config', 'set', 'backend', 'slurm'])
    res = cli_runner.invoke(main, ['config', 'get', 'backend'])
    assert res.output == 'slurm\n'

    cli_runner.invoke(main, ['config', 'unset', 'backend'])
    res = cli_runner.invoke(main, ['config', 'get', 'backend'])
    assert res.output == '\n'


def test_humanbool():
    assert humanbool('yes')
    assert humanbool('true')
    assert not humanbool('no')
    assert not humanbool('false')

    with pytest.raises(TypeError):
        humanbool('foo')


def test_cast_value():
    assert cast_value('yes')
    assert cast_value('true')
    assert not cast_value('no')
    assert not cast_value('false')
    assert cast_value('10') == 10
    assert cast_value('foo') == 'foo'

