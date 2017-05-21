from unittest import TestCase

from gwf.conf import Config, FileConfig


class TestConfig(TestCase):

    def test_length(self):
        c1 = Config(initial={'foo': 'bar'})
        self.assertEqual(len(c1), 1)

        c2 = Config(initial={'foo': 'bar', 'bar': 'baz'})
        self.assertEqual(len(c2), 2)

    def test_get_existing_key(self):
        c1 = Config(initial={'foo': 'bar', 'bar': 'baz'})
        self.assertEqual(c1['foo'], 'bar')
        self.assertEqual(c1['bar'], 'baz')

        c2 = Config(initial={'bar': {'foo': 'baz'}})
        self.assertEqual(c2['bar.foo'], 'baz')

        c3 = Config(initial={'bar': {'foo': {'foo': 'baz'}}})
        self.assertEqual(c3['bar.foo.foo'], 'baz')

    def test_get_nonexisting_key(self):
        c1 = Config(initial={'foo': {'baz': 'bar'}})
        with self.assertRaises(KeyError):
            c1['bar']

        with self.assertRaises(KeyError):
            c1['foo.bar']

        with self.assertRaises(KeyError):
            c1['foo.baz.qoox']

    def test_get_existing_key_with_default(self):
        c1 = Config(initial={'foo': {'baz': 'bar'}})
        self.assertEqual(c1.get('foo', 'baz'), {'baz': 'bar'})
        self.assertEqual(c1.get('foo.baz', 'baz'), 'bar')

    def test_get_nonexisting_key_with_default(self):
        c1 = Config(initial={'foo': {'baz': 'bar'}})
        self.assertEqual(c1.get('bar', 'baz'), 'baz')
        self.assertEqual(c1.get('bar', 'foo'), 'foo')

    def test_set_key(self):
        c1 = Config(initial={'foo': {'baz': 'bar'}})

        c1['bar'] = 'baz'
        self.assertEqual(c1['bar'], 'baz')

        c1['foo.bar'] = 'baz'
        self.assertEqual(c1['foo.bar'], 'baz')

        c1['foo.quux.foo.baz'] = 'tada'
        self.assertEqual(c1['foo.quux.foo.baz'], 'tada')

    def test_del_key(self):
        c1 = Config(initial={'foo': {'baz': 'bar'}})
        del c1['foo.baz']
        self.assertEqual(len(c1), 1)
        with self.assertRaises(KeyError):
            c1['foo.baz']

        del c1['foo']
        self.assertEqual(len(c1), 0)
        with self.assertRaises(KeyError):
            c1['foo']

    def test_iterkeys_1(self):
        c1 = Config(initial={'foo': 'bar', 'baz': 'foo'})
        self.assertEqual(set(c1.iterkeys()), {'foo', 'baz'})

    def test_iterkeys_2(self):
        c1 = Config(initial={'foo': {'bar': 'baz'}, 'baz': 'foo'})
        self.assertEqual(set(c1.iterkeys()), {'foo.bar', 'baz'})

    def test_iter(self):
        c1 = Config(initial={'foo': 'bar', 'baz': 'foo'})
        self.assertEqual(set(c1), {'foo', 'baz'})


def test_load_config_from_file(tmpdir):
    config_path = tmpdir.join('.gwfconf.json')
    config_path.write('{"foo": "bar"}\n')

    c = FileConfig.load(str(config_path))
    assert len(c) == 1
    assert c['foo'] == 'bar'


def test_load_config_from_nonexisting_file(tmpdir):
    config_path = tmpdir.join('.gwfconf.json')
    c = FileConfig.load(str(config_path))
    assert len(c) == 0


def test_dump_config_to_file(tmpdir):
    config_path = tmpdir.join('.gwfconf.json')
    config_path.write('{"foo": "bar"}\n')

    c1 = FileConfig.load(str(config_path))
    c1['baz'] = 'foo'
    c1.dump()

    c2 = FileConfig.load(str(config_path))
    assert len(c2) == 2
    assert c2['baz'] == 'foo'
    assert c2['foo'] == 'bar'
