import unittest
from unittest.mock import patch

from gwf.exceptions import TargetDoesNotExistError
from gwf import Target
from gwf.utils import parse_path, cache, get_file_timestamp, match_targets


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


class TestGetFileTimestamp(unittest.TestCase):

    @patch('gwf.utils.os.path.getmtime', return_value=42)
    def test_returns_modified_time_of_file_if_it_exists(self, mock_getmtime):
        self.assertEqual(get_file_timestamp('/some/file'), 42)
        mock_getmtime.assert_called_once_with('/some/file')

    @patch('gwf.utils.os.path.getmtime', side_effect=OSError)
    def test_returns_none_if_file_does_not_exist(self, mock_getmtime):
        self.assertIsNone(get_file_timestamp('/some/file'))
        mock_getmtime.assert_called_once_with('/some/file')


class TestMatchTargets(unittest.TestCase):

    def setUp(self):
        self.foobar = Target('foobar', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.foofoo = Target('foobar', inputs=[], outputs=[], options={}, working_dir='/some/dir')
        self.barbar = Target('foobar', inputs=[], outputs=[], options={}, working_dir='/some/dir')

        self.targets = {
            'foobar': self.foobar,
            'foofoo': self.foofoo,
            'barbar': self.barbar,
        }

    def test_match_one_target_by_name(self):
        self.assertEqual(match_targets(['foobar'], self.targets), {self.foobar})

    def test_match_two_targets_by_name(self):
        self.assertEqual(match_targets(['foobar', 'foofoo'], self.targets), {self.foobar, self.foofoo})

    def test_match_two_targets_by_pattern(self):
        self.assertEqual(match_targets(['foo*'], self.targets), {self.foobar, self.foofoo})

    def test_match_one_target_by_pattern(self):
        self.assertEqual(match_targets(['bar*'], self.targets), {self.barbar})

    def test_match_two_targets_by_pattern_and_one_by_name(self):
        self.assertEqual(match_targets(['foo*', 'barbar'], self.targets), {self.foobar, self.foofoo, self.barbar})

    def test_match_zero_targets_by_pattern(self):
        self.assertEqual(match_targets(['baz*'], self.targets), set())

    def test_raise_exception_if_target_name_does_not_exist(self):
        with self.assertRaises(TargetDoesNotExistError):
            match_targets(['baz'], self.targets)
