from __future__ import absolute_import, print_function

import imp
import inspect
import os.path
import sys

from .exceptions import GWFException


_target_repr = (
    '{}(name={!r}, inputs={!r}, outputs={!r}, options={!r}, working_dir={!r}, '
    'spec={!r})'
)


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


def _norm_paths(working_dir, paths):
    return [_norm_path(working_dir, path) for path in paths]


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

    def __init__(self):
        self.targets = {}

        filename = inspect.getfile(sys._getframe(1))
        self.working_dir = os.path.dirname(os.path.realpath(filename))

    def _add_target(self, target):
        if target.name in self.targets:
            raise GWFException(
                'Target "{}" already exists in workflow.'.format(target.name)
            )

        self.targets[target.name] = target

    def target(self, name, inputs, outputs, **options):
        """Create a target."""
        new_target = Target(
            name, inputs, outputs, options, working_dir=self.working_dir
        )

        self._add_target(new_target)
        return new_target

    def include_path(self, path):
        """Include targets of another workflow into this workflow."""
        other_workflow = _import_object(path)
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
