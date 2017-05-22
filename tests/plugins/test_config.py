from gwf.cli import main
from tests import CliTestCase

from gwf.plugins.config import humanbool, cast_value


class TestConfig(CliTestCase):

    def test_set_get_unset(self):
        self.runner.invoke(main, ['config', 'set', 'backend', 'slurm'])
        res = self.runner.invoke(main, ['config', 'get', 'backend'])
        self.assertEqual(res.output, 'slurm\n')

        self.runner.invoke(main, ['config', 'unset', 'backend'])
        res = self.runner.invoke(main, ['config', 'get', 'backend'])
        self.assertEqual(res.output, '\n')

    def test_humanbool(self):
        self.assertTrue(humanbool('yes'))
        self.assertTrue(humanbool('true'))
        self.assertFalse(humanbool('no'))
        self.assertFalse(humanbool('false'))

        with self.assertRaises(TypeError):
            humanbool('foo')

    def test_cast_value(self):
        self.assertTrue(cast_value('yes'))
        self.assertTrue(cast_value('true'))
        self.assertFalse(cast_value('no'))
        self.assertFalse(cast_value('false'))
        self.assertEqual(cast_value('10'), 10)
        self.assertEqual(cast_value('foo'), 'foo')

