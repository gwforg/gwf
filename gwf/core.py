from __future__ import absolute_import, print_function

import inspect
import itertools
import logging
import os.path
import sys
from collections import defaultdict

from .exceptions import (CircularDependencyError,
                         FileProvidedByMultipleTargetsError,
                         FileRequiredButNotProvidedError, TargetExistsError,
                         IncludeWorkflowError)
from .utils import (cache, get_file_timestamp, import_object, iter_inputs,
                    iter_outputs, timer)

logger = logging.getLogger(__name__)


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
    """Represents a target."""

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
        format_str = (
            '{}(name={!r}, inputs={!r}, outputs={!r}, options={!r}, '
            'working_dir={!r}, spec={!r})'
        )

        return format_str.format(
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

    def __init__(self, name=None, working_dir=None):
        self.name = name
        self.targets = {}

        self.working_dir = working_dir
        if self.working_dir is None:
            # Get the frame object of whatever called the Workflow.__init__
            # and extract the path of the file which is was defined in. Then
            # normalize the path and get the directory of the file.
            filename = inspect.getfile(sys._getframe(1))
            self.working_dir = os.path.dirname(os.path.realpath(filename))

    def _add_target(self, target):
        if target.name in self.targets:
            raise TargetExistsError(target)

        self.targets[target.name] = target

    def target(self, name, inputs=None, outputs=None, **options):
        """Create a target and add it to the :class:`gwf.core.Workflow`."""

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
        return other_workflow

    def include_workflow(self, other_workflow, namespace=None):
        """Include targets from another `Workflow` object in this workflow."""
        if other_workflow.name is None and namespace is None:
            raise IncludeWorkflowError(
                'The included workflow has not been assigned a name. To '
                'include the workflow you must set the namespace argument.'
            )
        namespace_prefix = namespace or other_workflow.name
        if namespace_prefix == self.name:
            raise IncludeWorkflowError(
                'The included workflow has the same name as this workflow'
            )

        for target in other_workflow.targets.values():
            target.name = '{}.{}'.format(namespace_prefix, target.name)
            self._add_target(target)

    def include(self, other_workflow, namespace=None):
        """Include another workflow into this workflow.

        This method can be given either an :class:`gwf.core.Workflow` instance,
        a module or a path to a workflow file.

        If a module or path the workflow object to include will be determined
        according to the following rules:

        1. If a module object is given, the module must define an attribute
           named `gwf` containing a :class:`gwf.core.Workflow` object.
        2. If a path is given it must point to a file defining a module with an
           attribute named `gwf` containing a :class:`gwf.core.Workflow`
           object.

        When a :class:`gwf.core.Workflow` instance has been obtained, all
        targets will be included directly into this workflow. To avoid name
        clashes the `namespace` argument must be provided. For example::

            workflow1 = Workflow()
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2, namespace='wf1')

        The workflow now contains two targets named `TestTarget` (defined in
        `workflow2`) and `wf1.TestTarget` (defined in `workflow1`). The
        `namespace` parameter can be left out if the workflow to be included
        has been named::

            workflow1 = Workflow(name='wf1')
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2)

        This yields the same result as before. The `namespace` argument can be
        used to override the specified name::

            workflow1 = Workflow(name='wf1')
            workflow1.target('TestTarget')

            workflow2 = Workflow()
            workflow2.target('TestTarget')

            workflow1.include(workflow2, namespace='foo')

        The workflow will now contain targets named `TestTarget` and
        `foo.TestTarget`.
        """
        if isinstance(other_workflow, Workflow):
            self.include_workflow(other_workflow)
        elif isinstance(other_workflow, str):
            self.include_path(other_workflow)
        elif inspect.ismodule(other_workflow):
            self.include_workflow(getattr(other_workflow, 'gwf'))
        else:
            raise TypeError('First argument must be either a string or a '
                            'Workflow object.')

    def __repr__(self):
        return '{}(name={!r}, working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.name, self.working_dir, self.targets
        )


class PreparedWorkflow:

    """Represents a finalized workflow graph."""

    def __init__(self, workflow=None):
        self.targets = {}
        self.working_dir = None

        if workflow is not None:
            self.prepare(workflow)

    def prepare(self, workflow):
        """Prepare this workflow given a :class:`gwf.core.Workflow` instance."""
        self.targets = workflow.targets
        self.working_dir = workflow.working_dir

        logger.debug(
            'preparing workflow with %s targets defined.',
            len(self.targets)
        )

        self.provides = self.prepare_file_providers()
        self.dependencies = self.prepare_dependencies(self.provides)
        self.dependents = self.prepare_dependents(self.dependencies)
        self._check_for_circular_dependencies()

        self.file_cache = self.prepare_file_cache()

    @timer('prepared file providers in %.3fms.', logger=logger)
    def prepare_file_providers(self):
        provides = {}
        for target, path in iter_outputs(self.targets.values()):
            if path in provides:
                raise FileProvidedByMultipleTargetsError(
                    path, provides[path].name, target
                )

            provides[path] = target
        return provides

    @timer('prepared dependencies in %.3fms.', logger=logger)
    def prepare_dependencies(self, rovides):
        dependencies = defaultdict(list)
        for target, path in iter_inputs(self.targets.values()):
            if os.path.exists(path):
                continue

            if path not in self.provides:
                raise FileRequiredButNotProvidedError(path, target)
            dependencies[target].append(self.provides[path])
        return dependencies

    @timer('prepared dependents in %.3fms.', logger=logger)
    def prepare_dependents(self, dependencies):
        dependents = defaultdict(list)
        for target, deps in self.dependencies.items():
            for dep in deps:
                dependents[dep].append(target)
        return dependents

    @timer('prepared file cache in %.3fms.', logger=logger)
    def prepare_file_cache(self):
        input_iter = iter_inputs(self.targets.values())
        output_iter = iter_outputs(self.targets.values())
        return {path: get_file_timestamp(path)
                for _, path in itertools.chain(input_iter, output_iter)}

    @timer('checked for circular dependencies in %.3fms.', logger=logger)
    def _check_for_circular_dependencies(self):
        for target in self.targets.values():
            for dep in self.dependencies[target]:
                if target in _get_deep_dependencies(dep, self.dependencies):
                    raise CircularDependencyError(target)

    @cache
    def should_run(self, target):
        if any(self.should_run(dep) for dep in self.dependencies[target]):
            logger.debug(
                '%s should run because one of its dependencies should run.',
                target.name
            )
            return True

        if target.is_sink:
            logger.debug(
                '%s should run because it is a sink.',
                target.name
            )
            return True

        if any(self.file_cache[path] is None for path in target.outputs):
            logger.debug(
                '%s should run because one of its output files does not exist.',
                target.name
            )
            return True

        if target.is_source:
            logger.debug(
                '%s should not run because it is a source.',
                target.name
            )
            return False

        youngest_in_ts, youngest_in_path = max(
            self.file_cache[path] for path in target.inputs
        )

        logger.debug(
            '%s is the youngest input file of %s with timestamp %s.',
            youngest_in_path,
            target.name,
            youngest_in_ts
        )

        oldest_out_ts, oldest_out_path = min(
            self.file_cache[path] for path in target.outputs
        )

        logger.debug(
            '%s is the oldest output file of %s with timestamp %s.',
            oldest_out_path,
            target.name,
            youngest_in_ts
        )

        logger.debug(
            '%s should run since %s is larger than %s.',
            target.name,
            youngest_in_ts,
            oldest_out_ts
        )

        return youngest_in_ts > oldest_out_ts

    def __repr__(self):
        return '{}(working_dir={!r}, targets={!r})'.format(
            self.__class__.__name__, self.working_dir, self.targets
        )
