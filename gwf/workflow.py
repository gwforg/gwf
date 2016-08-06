from __future__ import absolute_import, print_function

import imp
import inspect
import os.path
import sys
from collections import defaultdict

from . import GWFException
from .utils import cache


ex_msg_file_provided_by_multple_targets = '''
    File "{}" provided by multiple targets "{}" and "{}".
'''.strip()

ex_msg_file_required_but_not_provided = '''
    File "{}" is required by "{}", but does not exist and is not provided by \
    a target.
'''.strip()

target_repr = '''
    {}(name={!r}, inputs={!r}, outputs={!r}, options={!r}, working_dir={!r}, spec={!r})
'''.strip()


def _import_object(path, default_obj='gwf'):
    if not os.path.isabs(path):
        path = os.path.abspath(os.path.join(os.getcwd(), path))

    comps = path.rsplit(':')
    if len(comps) == 2:
        path, obj = comps
    elif len(comps) == 1:
        path, obj = comps[0], default_obj
    else:
        raise ValueError('Invalid path.')

    basedir, filename = os.path.split(path)
    filename, ext = os.path.splitext(filename)

    mod_loc = imp.find_module(filename, [basedir])
    mod = imp.load_module(filename, *mod_loc)
    return getattr(mod, obj)


def _norm_path(working_dir, path):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


@cache
def _get_file_timestamp(filename):
    return os.path.getmtime(filename)


def dependencies(nodes, target_name):
    """Return all tasks necessary for building the target.

    The set of tasks is just returned as set.
    """
    root = nodes[target_name]

    # Working with a list to preserve the order. It makes lookups slower but
    # hopefully these sets won't be terribly long ... if it becomes a problem
    # it is easy enough to fix it.
    processed = []

    def dfs(node):
        if node in processed:
            return
        else:
            for dep in node.depends_on:
                dfs(dep)
            processed.append(node)

    dfs(root)
    return processed


def get_execution_schedule(nodes, target_name):
    """Linearize the targets to be run.

    Returns a list of tasks to be run (in the order they should run or
    be submitted to the cluster to make sure dependencies are handled
    correctly) and a set of the names of tasks that will be scheduled
    (to make sure dependency flags are set in the submission command).
    """

    root = nodes[target_name]

    # If the target is already in the queue we just dismiss the scheduling
    # right away... this because we need to handle dependent nodes in the
    # queue differently, since for those we need wait for completion.
    if root.job_in_queue:
        return [], set()

    processed = set()
    scheduled = set()
    job_schedule = []

    def dfs(node):
        if node in processed:
            # we have already processed the node, and
            # if we should run the target name is scheduled
            # otherwise it isn't.
            return node.target.name in scheduled

        # schedule all dependencies before we schedule this task
        for dep in node.depends_on:
            dfs(dep)

        # If this task needs to run, then schedule it
        if node.job_in_queue or node.should_run:
            job_schedule.append(node)
            scheduled.add(node.name)

        processed.add(node)

    dfs(root)

    return job_schedule, scheduled


class Target(object):

    def __init__(self, name, inputs, outputs, options, working_dir, spec=None):
        self.name = name

        self.inputs = [_norm_path(working_dir, path) for path in inputs]
        self.outputs = [_norm_path(working_dir, path) for path in outputs]

        self.options = options
        self.working_dir = working_dir

        self.spec = spec

    @property
    def is_source(self):
        """Return whether this target is a source.

        A target is a source if it does not depend on any files."""
        return not self.inputs

    @property
    def is_sink(self):
        """Return whether this target is a sink.

        A target is a sink if it does not output any files.
        """
        return not self.outputs

    def __lshift__(self, spec):
        self.spec = spec

    def __repr__(self):
        return target_repr.format(
            self.__class__.__name__,
            self.name,
            self.inputs,
            self.outputs,
            self.options,
            self.working_dir,
            self.spec,
        )

    def __str__(self):
        return self.name


class Workflow(object):

    def __init__(self):
        self._backend = None
        self._targets = {}
        self._dependency_graph = {}

        filename = inspect.getfile(sys._getframe(1))
        self.working_dir = os.path.dirname(os.path.realpath(filename))

    def _add_target(self, target):
        if target.name in self._targets:
            raise GWFException(
                'Target "{}" already exists in workflow.'.format(target.name)
            )

        self._targets[target.name] = target

    def target(self, name, inputs, outputs, **options):
        """Create a target."""
        new_target = Target(name, inputs, outputs, options,
                            working_dir=self.working_dir)

        self._add_target(new_target)
        return new_target

    def include_path(self, path):
        """Include targets of another workflow into this workflow."""
        other_workflow = _import_object(path)
        self.include_workflow(other_workflow)

    def include_workflow(self, other_workflow):
        """Include targets from another `Workflow` object in this workflow."""
        for target in other_workflow._targets.values():
            self._add_target(target)

    def include(self, path_or_workflow):
        """Include another workflow into this workflow."""
        if isinstance(path_or_workflow, Workflow):
            self.include_workflow(path_or_workflow)
        elif isinstance(path_or_workflow, str):
            self.include_path(path_or_workflow)
        elif inspect.ismodule(path_or_workflow):
            self.include_workflow(getattr(path_or_workflow, 'gwf'))
        else:
            raise TypeError('First argument must be either a string or a '
                            'Workflow object.')

    def prepare(self):
        """Prepare workflow to run.

        After the workflow has been prepared, three dictionaries will be
        available on the object:

        `provides`:
            a dictionary of file paths and the corresponding targets that
            produce the files.

        `dependencies`:
            a dictionary where each key is a `Target` and the value is a list
            of `Target`s that the key `Target` depends on.

        `dependents`:
            a dictionary where each key is a `Target` and the value is a list
            of `Targets`s that depend on the key `Target`.
        """
        provides = {}
        dependencies = defaultdict(list)
        dependents = defaultdict(list)

        for target in self._targets.values():
            for path in target.outputs:
                if path in provides:
                    raise GWFException(
                        ex_msg_file_provided_by_multple_targets.format(
                            path, provides[path].name, target.name
                        )
                    )

                provides[path] = target

        for target in self._targets.values():
            for path in target.inputs:
                if os.path.exists(path):
                    continue

                if path not in provides:
                    raise GWFException(
                        ex_msg_file_required_but_not_provided.format(
                            path, target.name
                        )
                    )

                dependencies[target].append(provides[path])

        for target, deps in dependencies.items():
            for dep in deps:
                dependents[dep].append(target)

        self.provides = provides
        self.dependencies = dependencies
        self.dependents = dependents

    @cache
    def should_run(self, target):
        if any(self.should_run(dep) for dep in self.dependencies[target]):
            return True

        if target.is_sink:
            return True

        if any(not os.path.exists(path) for path in target.outputs):
            return True

        if target.is_source:
            return False

        youngest_in_ts, youngest_in_path = \
            max((_get_file_timestamp(path), path) for path in target.inputs)

        oldest_out_ts, oldest_out_path = \
            min((_get_file_timestamp(path), path) for path in target.outputs)

        return youngest_in_ts > oldest_out_ts

    def run(self, backend):
        pass
