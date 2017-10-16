import unittest

from gwf.utils import parse_path, cache


class TestCache(unittest.TestCase):

    def test_returns_same_object_when_called_twice_with_same_args(self):
        func = cache(lambda x: object())
        obj1 = func(42)
        obj2 = func(42)
        self.assertEqual(id(obj1), id(obj2))

    def test_does_not_return_same_object_when_called_with_diff_args(self):
        func = cache(lambda x: object())
        obj1 = func(42)
        obj2 = func(43)
        self.assertNotEqual(id(obj1), id(obj2))


class TestParsePath(unittest.TestCase):

    def test_without_object_specified(self):
        basedir, filename, obj = parse_path('/some/dir/workflow.py', 'workflow_obj')
        self.assertEqual(filename, 'workflow')
        self.assertEqual(basedir, '/some/dir')
        self.assertEqual(obj, 'workflow_obj')

    def test_with_object_specified(self):
        basedir, filename, obj = parse_path('/some/dir/workflow.py:other_obj', 'workflow_obj')
        self.assertEqual(filename, 'workflow')
        self.assertEqual(basedir, '/some/dir')
        self.assertEqual(obj, 'other_obj')

    def test_parse_invalid_path_raises_value_error(self):
        with self.assertRaises(ValueError):
            parse_path('/some/dir/workflow.py::other_obj', 'workflow_obj')
