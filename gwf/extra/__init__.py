import glob as _glob
import inspect
import marshal
import os
import os.path
import shelve
import subprocess
import sys

import gwf

from gwf.helpers import file_exists, get_file_timestamp, make_absolute_path

def _list(x):
    """Wrap x as a singleton in a list if it isn't a list already."""
    if hasattr(x, '__iter__'):
        return list(x)
    else:
        return [x]


## Useful helper functions...
def glob(pattern):
    """Returns a list of filenames matching the shell glob `pattern'."""
    if pattern.startswith('/'):
        glob_pattern = pattern
    else:
        filename = inspect.getfile(sys._getframe(1))
        working_dir = os.path.dirname(os.path.realpath(filename))
        glob_pattern = os.path.join(working_dir, pattern)
    return _glob.glob(glob_pattern)


def shell(command):
    """Execute the shell command `command' and return the resulting output as a list of tokens."""
    return subprocess.check_output(command, shell=True).split()


## sub-workflows
def include_workflow(path):
    """Import targets from another workflow file.

    If called with a path that ends in ".py" it is assumed to be the name of
    a workflow file and that file is included. If the path doesn't end in ".py"
    it is assumed to be a directory and "/workflow.py" is appended to it.
    """

    if path.endswith(".py"):
        workflow_file = path
    else:
        workflow_file = path + '/workflow.py'

    if not workflow_file.startswith('/'):
        # include relative to the workflow that includes...
        filename = inspect.getfile(sys._getframe(1))
        working_dir = os.path.dirname(os.path.realpath(filename))
        workflow_file = os.path.join(working_dir, workflow_file)

    execfile(workflow_file)


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
            options = dict(spec[0].items() + options.items())
            spec = spec[1]

        if self.name in gwf.ALL_TARGETS:
            print 'Warning: Target', self.name, 'defined more than once.'

        new_target = gwf.Target(self.name, options, spec)
        gwf.ALL_TARGETS[self.name] = new_target


class _memorize_wrapper(object):
    def __init__(self, func, options):
        self.func = func
        self.options = options
        if 'input' not in self.options:
            self.options['input'] = []
        self.options['input'] = _list(self.options['input'])

        filename = inspect.getfile(sys._getframe(2))
        self.working_dir = os.path.dirname(os.path.realpath(filename))

        # The database of remembered results
        memory_dir = self.memory_dir()
        db_file = '{}/{}'.format(memory_dir, self.func.func_name)
        self.results_db = shelve.open(db_file, writeback=True)
        flag_file = '{}/{}.flag'.format(memory_dir, self.func.func_name)

        if self.should_run(flag_file):
            # Clear the database of existing results since these must be out of date
            if os.path.exists(flag_file):
                os.unlink(flag_file)

            self.results_db.clear()
            open(flag_file, 'w').close() # touch the flag

        # Remember the byte code for the function so it gets run again if it changes
        current_code_string = marshal.dumps(func.func_code)
        if '-code-' not in self.results_db or current_code_string != self.results_db['-code-']:
            # nuke the old and create a new ... the old base based on old code
            if os.path.exists(flag_file):
                os.unlink(flag_file)
            self.results_db.clear()
            self.results_db = shelve.open(db_file)
            self.results_db['-code-'] = current_code_string
            open(flag_file, 'w').close() # touch the flag

    def memory_dir(self):
        memory_directory = os.path.join(self.working_dir, '.memory')
        if not os.path.exists(memory_directory):
            os.makedirs(memory_directory)
        return memory_directory

    def __call__(self, *args):
        if not str(args) in self.results_db:
            self.results_db[str(args)] = self.func(*args)
        return self.results_db[str(args)]

    def should_run(self, output_file):
        """Test if this target needs to be run based on whether input
        and output files exist and on their time stamps. Doesn't check
        if upstream targets need to run, only this task; upstream tasks
        are handled by the dependency graph. """

        if not file_exists(make_absolute_path(self.working_dir, output_file)):
            return True

        for infile in self.options['input']:
            if not file_exists(make_absolute_path(self.working_dir, infile)):
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
            timestamp = get_file_timestamp(make_absolute_path(self.working_dir, infile))
            if youngest_in_timestamp is None \
                    or youngest_in_timestamp < timestamp:
                youngest_in_timestamp = timestamp
        assert youngest_in_timestamp is not None

        out_timestamp = get_file_timestamp(make_absolute_path(self.working_dir, output_file))
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




_FUNCTION_TEMPLATE_TEMPLATE = """

import marshal
import shelve
import types
import sys


with open("{marshal_file}") as marshal_file:
    code = marshal.loads(marshal_file.read())
func = types.FunctionType(code, globals(), "{function_name}")

args_db = shelve.open("{arguments_db_file}")
args = args_db[sys.argv[1]]

func(*args)

"""

class _function_template_wrapper(object):

    def marshal_dir(self):
        d = os.path.join(self.working_dir, '.marshal/')
        if not os.path.exists(d):
            os.makedirs(d)
        return d

    def __init__(self, func, options):
        self.func = func
        self.options = options

        self.working_dir = os.path.dirname(os.path.realpath(inspect.getfile(sys._getframe(2))))

        self.marshal_file = os.path.join(self.marshal_dir(), '{}.code'.format(func.func_name))
        self.arguments_db_file = os.path.join(self.marshal_dir(), '{}.args'.format(func.func_name))
        self.python_file = os.path.join(self.marshal_dir(), '{}.py'.format(func.func_name))

        with open(self.marshal_file, 'w') as marshal_file:
            marshal_file.write(marshal.dumps(func.func_code))
        self.arguments_db = shelve.open(self.arguments_db_file)
        with open(self.python_file, 'w') as python_file:
            script = _FUNCTION_TEMPLATE_TEMPLATE.format(marshal_file=self.marshal_file,
                                                        function_name=self.func.func_name,
                                                        arguments_db_file=self.arguments_db_file,
                                                        )
            python_file.write(script)

        self.args_counter = 0

    def __call__(self, *args):

        self.arguments_db["{}".format(self.args_counter)] = args
        self.args_counter += 1

        options = dict(self.options.items())
        spec = '''
        python {python_file} {argument_index}
        '''.format(python_file=self.python_file, argument_index=self.args_counter-1)

        return options, spec

class function_template(object):

    def __init__(self, **options):
        self.options = options

    def __call__(self, func):
        return _function_template_wrapper(func, self.options)


class _Infix:
    def __init__(self, function):
        self.function = function
    def __rlshift__(self, other):
        return _Infix(lambda x, self=self, other=other: self.function(other, x))
    def __rshift__(self, other):
        return self.function(other)

def splice_path(path, directory=None, suffix=None, tag=None):
    orig_dir, base_name = os.path.split(path)
    base_name, orig_suffix = os.path.splitext(base_name)
    if not directory:
        directory = orig_dir
    if not suffix:
        suffix = orig_suffix
    if tag:
        base_name = "{}_{}".format(tag, base_name)
    return os.path.join(directory, base_name) + suffix

tag = _Infix(lambda l, tag: [splice_path(x, tag=tag) for x in l])
suffix = _Infix(lambda l, suffix: [splice_path(x, suffix=suffix) for x in l])
outdir = _Infix(lambda l, outdir: [splice_path(x, directory=outdir) for x in l])
