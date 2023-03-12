import logging
import os.path

import attrs

from .core import CachedFilesystem, Target, hash_spec
from .exceptions import WorkflowError

logger = logging.getLogger(__name__)


class Reason:
    @attrs.frozen
    class _Base:
        def reason(self):
            return NotImplementedError("reason()")

        def __str__(self):
            return self.reason()

    @attrs.frozen
    class IsSink(_Base):
        target: Target
        scheduled: bool = True

        def reason(self):
            return "is a sink"

        def long_reason(self):
            return self.reason()

    @attrs.frozen
    class IsSource(_Base):
        target: Target
        scheduled: bool = False

        def reason(self):
            return "not scheduled because it is a source"

        def long_reason(self):
            return self.reason()

    @attrs.frozen
    class DependencyScheduled(_Base):
        target: Target = attrs.field()
        dependencies: list = attrs.field()
        scheduled: bool = True

        def reason(self):
            return "a dependency was scheduled"

        def long_reason(self):
            return (
                f"scheduled because the dependency {self.dependencies[0]} was scheduled"
            )

    @attrs.frozen
    class MissingOutput(_Base):
        target: Target = attrs.field()
        path: str = attrs.field()
        scheduled: bool = True

        def reason(self):
            return "an output file is missing"

        def long_reason(self):
            path = os.path.relpath(self.path)
            return f"output file {path} is missing"

    @attrs.frozen
    class OutOfDate(_Base):
        target: Target = attrs.field()
        input_path: str = attrs.field()
        output_path: str = attrs.field()
        scheduled: bool = True

        def reason(self):
            return "is out-of-date"

        def long_reason(self):
            input_path = os.path.relpath(self.input_path)
            output_path = os.path.relpath(self.output_path)
            return (
                f"scheduled because input file {input_path} is newer "
                f"than output file {output_path}"
            )

    @attrs.frozen
    class SpecChanged(_Base):
        target: Target = attrs.field()
        old_hash: str = attrs.field()
        curr_hash: str = attrs.field()
        scheduled: bool = True

        def reason(self):
            return "spec has changed"

        def long_reason(self):
            return self.reason()

    @attrs.frozen
    class UpToDate(_Base):
        target: Target = attrs.field()
        scheduled: bool = False

        def reason(self):
            return "is up-to-date"

        def long_reason(self):
            return self.reason()


def linearize_plan(plan):
    linear_plan = []
    visited = set()

    def traverse(reason):
        deps = list()
        if reason.scheduled and isinstance(reason, Reason.DependencyScheduled):
            for dep in reason.dependencies:
                if dep.scheduled:
                    deps.append(dep.target)
                traverse(dep)

        if reason.target.name not in visited:
            linear_plan.append((reason, deps))
            visited.add(reason.target.name)

    traverse(plan)
    return linear_plan


def schedule_workflow(graph, fs=None, spec_hashes=None, endpoints=None):
    if endpoints is None:
        endpoints = graph.endpoints()

    if fs is None:
        fs = CachedFilesystem()

    reasons = dict()

    def get_reason(target, graph, fs, spec_hashes):
        scheduled_deps = []
        for dep in graph.dependencies[target]:
            reason = schedule_cached(dep)
            if reason.scheduled:
                scheduled_deps.append(reason)

        if scheduled_deps:
            return Reason.DependencyScheduled(target, scheduled_deps)

        if spec_hashes is not None:
            curr_hash = hash_spec(target.spec)
            old_hash = spec_hashes.get(target.name)
            if curr_hash != old_hash:
                return Reason.SpecChanged(target, old_hash, curr_hash)

        # Check whether all input files actually exists are are being provided
        # by another target. If not, it's an error.
        for path in target.flattened_inputs():
            if path in graph.unresolved and not fs.exists(path):
                msg = (
                    'File "{}" is required by "{}", but does not exist and is not '
                    "provided by any target in the workflow."
                ).format(path, target)
                raise WorkflowError(msg)

        for path in target.flattened_outputs():
            if not fs.exists(path):
                return Reason.MissingOutput(target, path)

        if not target.inputs:
            return Reason.IsSource(target)

        youngest_in_ts, youngest_in_path = max(
            (fs.changed_at(path), path) for path in target.flattened_inputs()
        )
        logger.debug(
            "%s is the youngest input file of %s with timestamp %s",
            youngest_in_path,
            target,
            youngest_in_ts,
        )

        if not target.outputs:
            return Reason.IsSink(target)

        oldest_out_ts, oldest_out_path = min(
            (fs.changed_at(path), path) for path in target.flattened_outputs()
        )
        logger.debug(
            "%s is the oldest output file of %s with timestamp %s",
            oldest_out_path,
            target,
            oldest_out_ts,
        )

        if youngest_in_ts > oldest_out_ts:
            return Reason.OutOfDate(target, youngest_in_path, oldest_out_path)
        return Reason.UpToDate(target)

    def schedule_cached(target):
        if target not in reasons:
            reason = get_reason(target, graph, fs, spec_hashes)
            reasons[target] = reason
        return reasons[target]

    for endpoint in endpoints:
        schedule_cached(endpoint)

    return reasons
