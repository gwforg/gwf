from __future__ import absolute_import, print_function

import inspect
import os.path
import sys
from collections import defaultdict

from .exceptions import (CircularDependencyError,
                         FileProvidedByMultipleTargetsError,
                         FileRequiredButNotProvidedError, TargetExistsError)
from .utils import import_object, iter_inputs, iter_outputs, get_file_timestamp

_target_repr = (
    '{}(name={!r}, inputs={!r}, outputs={!r}, options={!r}, working_dir={!r}, '
    'spec={!r})'
)


def _norm_path(working_dir, path):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(working_dir, path))


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


def _get_deep_dependencies(target, dependencies):
    """Return all `Target`s necessary for building `target`.

    The set of tasks is just returned as set.
    """
    # Working with a list to preserve the order. It makes lookups slower
    # but hopefully these sets won't be terribly long ... if it becomes a
    # problem it is easy enough to fix it.
    processed = []

    def dfs(other_target):
        if other_target in processed:
            return
        else:
            processed.append(other_target)
            for dep in dependencies[other_target]:
                dfs(dep)

    dfs(target)
    return processed


class Target(object):

    def __init__(self, name, inputs, outputs, options, working_dir, spec=None):
        self.name = name

        self.inputs = _norm_paths(working_dir, inputs)
        self.outputs = _norm_paths(working_dir, outputs)

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
        return self

    def __repr__(self):
        return _target_repr.format(
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

    """Represents a workflow."""

    def __init__(self, working_dir=None):
        self.targets = {}

        self.working_dir = working_dir
        if self.working_dir is None:
            # Get the frame object of whatever called the Workflow.__init__
            # and extract the path of the file which is was defined in. Then
            # normalize the path and get the directory of the file.
            #
            # TODO: Figure out whether this can be replaced with a simple
            # os.getcwd() call.
            filename = inspect.getfile(sys._getframe(1))
            self.working_dir = os.path.dirname(os.path.realpath(filename))

    def _add_target(self, target):
        if target.name in self.targets:
            raise TargetExistsError(target)

        self.targets[target.name] = target

    def target(self, name, inputs=None, outputs=None, **options):
        """Create a target and add it to the `Workflow`."""

        if inputs is None:
            inputs = []
        if outputs is None:
            outputs = []

        new_target = Target(
            name, inputs, outputs, options, working_dir=self.working_dir
        )

        self._add_target(new_target)
        return new_target

    def include_path(self, path):
        """Include targets of another workflow into this workflow."""
        other_workflow = import_object(path)
        self.include_workflow(other_workflow)

    def include_workflow(self, other_workflow):
        """Include targets from another `Workflow` object in this workflow."""
        for target in other_workflow.targets.values():
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

    def __repr__(self):
        return '{}(working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.working_dir, self.targets
        )


class PreparedWorkflow:

    def __init__(self, workflow=None):
        self.targets = {}
        self.working_dir = None

        if workflow is not None:
            self.prepare(workflow)

    def prepare(self, workflow):
        self.targets = workflow.targets
        self.working_dir = workflow.working_dir

        self.provides = self.prepare_file_providers()
        self.dependencies = self.prepare_dependencies(self.provides)
        self.dependents = self.prepare_dependents(self.dependencies)

        self._check_for_circular_dependencies()

        self.file_cache = self.prepare_file_cache()

    def prepare_file_providers(self):
        provides = {}
        for target, path in iter_outputs(self.targets.values()):
            if path in provides:
                raise FileProvidedByMultipleTargetsError(
                    path, provides[path].name, target
                )

            provides[path] = target
        return provides

    def prepare_dependencies(self, rovides):
        dependencies = defaultdict(list)
        for target, path in iter_inputs(self.targets.values()):
            if os.path.exists(path):
                continue

            if path not in self.provides:
                raise FileRequiredButNotProvidedError(path, target)
            dependencies[target].append(self.provides[path])
        return dependencies

    def prepare_dependents(self, dependencies):
        dependents = defaultdict(list)
        for target, deps in self.dependencies.items():
            for dep in deps:
                dependents[dep].append(target)
        return dependents

    def _check_for_circular_dependencies(self):
        for target in self.targets.values():
            for dep in self.dependencies[target]:
                if target in _get_deep_dependencies(dep, self.dependencies):
                    raise CircularDependencyError(target)

    def __repr__(self):
        return '{}(working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.working_dir, self.targets
        )
