import unittest
from unittest.mock import patch

from gwf.exceptions import GWFError
from gwf.utils import (_split_import_path, cache, get_file_timestamp,
                       load_workflow)


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


class TestImportObject(unittest.TestCase):

    def test_split_import_path_without_object_specified(self):
        filename, basedir, obj = _split_import_path(
            '/some/dir/workflow.py', 'workflow_obj'
        )

        self.assertEqual(filename, 'workflow')
        self.assertEqual(basedir, '/some/dir')
        self.assertEqual(obj, 'workflow_obj')

    def test_split_import_path_with_object_specified(self):
        filename, basedir, obj = _split_import_path(
            '/some/dir/workflow.py:other_obj', 'workflow_obj'
        )

        self.assertEqual(filename, 'workflow')
        self.assertEqual(basedir, '/some/dir')
        self.assertEqual(obj, 'other_obj')

    def test_split_invalid_import_path_raises_value_error(self):
        with self.assertRaises(ValueError):
            _split_import_path('/some/dir/workflow.py::other_obj', 'workflow_obj')

    @patch('gwf.utils.os.path.exists', return_value=True)
    @patch('gwf.utils.os.getcwd', return_value='/some/dir')
    @patch('gwf.utils.imp.find_module', return_value=(None, '', ('', '', None)))
    @patch('gwf.utils.imp.load_module')
    @patch('gwf.utils._split_import_path', return_value=('/some/dir/this', 'workflow', 'gwf'))
    def test_import_with_non_absolute_path_normalizes_path_and_loads_module(
            self, mock_split_import_path, mock_load_module, mock_find_module, mock_getcwd, mock_exists):

        load_workflow('this/workflow.py')

        mock_split_import_path.assert_called_once_with(
            '/some/dir/this/workflow.py', 'gwf'
        )
        self.assertEqual(mock_find_module.call_count, 1)
        self.assertEqual(mock_load_module.call_count, 1)

    def test_trying_to_load_workflow_from_nonexisting_path_raises_exception(self):
        with self.assertRaisesRegex(GWFError, '.*this/workflow\.py.*'):
            load_workflow('this/workflow.py')


class TestGetFileTimestamp(unittest.TestCase):

    @patch('gwf.utils.os.path.getmtime', return_value=42)
    def test_returns_modified_time_of_file_if_it_exists(self, mock_getmtime):
        self.assertEqual(get_file_timestamp('/some/file'), 42)
        mock_getmtime.assert_called_once_with('/some/file')

    @patch('gwf.utils.os.path.getmtime', side_effect=OSError)
    def test_returns_none_if_file_does_not_exist(self, mock_getmtime):
        self.assertIsNone(get_file_timestamp('/some/file'))
        mock_getmtime.assert_called_once_with('/some/file')
