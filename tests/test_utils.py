import unittest

from gwf.utils import cache


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
