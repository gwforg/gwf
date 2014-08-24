import os
import os.path
import sys
import inspect
import glob as _glob
import subprocess
import cPickle

import gwf_workflow

## Useful helper functions...


def glob(pattern):
    if pattern.startswith('/'):
        glob_pattern = pattern
    else:
        filename = inspect.getfile(sys._getframe(1))
        working_dir = os.path.dirname(os.path.realpath(filename))
        glob_pattern = os.path.join(working_dir, pattern)
    return _glob.glob(glob_pattern)


def shell(command):
    return subprocess.check_output(command, shell=True).split()


## Templates
class template(object):
    def __init__(self, **options):
        self.spec = None
        self.options = options

    def __lshift__(self, spec):
        self.spec = spec
        return self

    def __call__(self, **substitutions):
        def substitute(s):
            if type(s) == str:
                return s.format(**substitutions)
            elif hasattr(s, '__iter__'):
                return [substitute(x) for x in s]
            else:
                return s

        formatted_options = [(key, substitute(val)) for key, val in self.options.items()]
        options = dict(formatted_options)
        return options, self.spec.format(**substitutions)


# Targets...
class target(object):
    def __init__(self, name, **options):
        self.name = name
        self.options = options

    def __lshift__(self, spec):
        options = self.options

        if isinstance(spec, tuple):
            options = dict(options.items() + spec[0].items())
            spec = spec[1]

        if self.name in gwf_workflow.ALL_TARGETS:
            print 'Warning: Target', self.name, 'defined more than once.'

        new_target = gwf_workflow.Target(self.name, options, spec)
        gwf_workflow.ALL_TARGETS[self.name] = new_target



# FIXME: these helper functions are copied. We should only have them written in one place

_remembered_files = {}  # cache to avoid too much stat'ing
_remembered_timestamps = {}  # Use this to avoid too many stats that slows down the script

def _file_exists(filename):
    if filename not in _remembered_files:
        _remembered_files[filename] = os.path.exists(filename)
    return _remembered_files[filename]


def _get_file_timestamp(filename):
    if filename not in _remembered_timestamps:
        _remembered_timestamps[filename] = os.path.getmtime(filename)
    return _remembered_timestamps[filename]


def _make_absolute_path(working_dir, filename):
    if os.path.isabs(filename):
        abspath = filename
    else:
        abspath = os.path.join(working_dir, filename)
    return os.path.normpath(abspath)

def _list(x):
    """Wrap x as a singleton in a list if it isn't a list already."""
    if hasattr(x, '__iter__'):
        return list(x)
    else:
        return [x]


class _memorize_wrapper(object):
    def __init__(self, func, options):
        self.func = func
        self.options = options
        if 'input' not in self.options:
            self.options['input'] = []
        self.options['input'] = _list(self.options['input'])

        filename = inspect.getfile(sys._getframe(2))
        self.working_dir = os.path.dirname(os.path.realpath(filename))

    def __call__(self, *args):
        output_file = '{}.result'.format(self.func.func_name)
        if self.should_run(output_file):
            self.result = self.func(*args)
            with open(output_file, 'w') as memory:
                cPickle.dump(self.result, memory)
        else:
            with open(output_file, 'r') as memory:
                self.result = cPickle.load(memory)
        return self.result

    def should_run(self, output_file):
        """Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. """

        if not _file_exists(_make_absolute_path(self.working_dir, output_file)):
            return True

        for infile in self.options['input']:
            if not _file_exists(_make_absolute_path(self.working_dir, infile)):
                print """
The memorized function `{func_name}' depends on an input file "{infile}"
that doesn't exist.

Memorized functions cannot depend on the output of targets since they
need to run at the time the workflow is evaluated.
                """.format(func_name=self.func.func_name, infile=infile)
                sys.exit(1)

        # If no file is missing, it comes down to the time stamps. If we
        # only have output and no input, we assume the output is up to
        # date. Touching files and adding input can fix this behaviour
        # from the user side but if we have a program that just creates
        # files we don't want to run it whenever someone needs that
        # output just because we don't have time stamped input.

        if len(self.options['input']) == 0:
            return False

        # if we have both input and output files, check time stamps

        youngest_in_timestamp = None
        for infile in self.options['input']:
            timestamp = _get_file_timestamp(_make_absolute_path(self.working_dir, infile))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        out_timestamp = _get_file_timestamp(_make_absolute_path(self.working_dir, output_file))
        # The youngest in should be older than the output
        if youngest_in_timestamp >= out_timestamp:
            return True
        else:
            return False


class memorize(object):
    """Remember the output of a function between evaluations of the workflow.

    Can be used together with `input' files so it is reevaluated if data the output depends on changes.
    """

    def __init__(self, **options):
        self.options = options

    def __call__(self, func):
        return _memorize_wrapper(func, self.options)
