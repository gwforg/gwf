import logging
from functools import partial, lru_cache

from .backends.base import BackendStatus
from .core import Status

logger = logging.getLogger(__name__)


SUBMITTED_STATES = (
    Status.SUBMITTED,
    Status.RUNNING,
    Status.SHOULDRUN,
    Status.FAILED,
    Status.CANCELLED,
)


def schedule(
    endpoints,
    graph,
    fs,
    spec_hashes,
    status_func,
    submit_func,
    force=False,
    no_deps=False,
):
    schedule_cache = {}

    def _is_outdated(input_paths, output_paths):
        youngest_in = max(
            (fs.changed_at(p) for p in input_paths if fs.exists(p)),
            default=float("-inf"),
        )
        oldest_out = min(
            (fs.changed_at(p) for p in output_paths),
            default=float("inf"),
        )
        return youngest_in > oldest_out

    @lru_cache(maxsize=None)
    def _all_downstream_outputs_exist(target):
        if not (dependents := graph.dependents[target]):
            return False
        for dep in dependents:
            if dep_non_temp := set(dep.flattened_outputs()) - graph.temporary:
                # Dependent has non-temp outputs – they must all exist and be
                # up-to-date for the path to be considered materialised.
                if not all(fs.exists(out) for out in dep_non_temp):
                    return False
                dep_non_temp_inputs = set(dep.flattened_inputs()) - graph.temporary
                if _is_outdated(dep_non_temp_inputs, dep_non_temp):
                    return False
            elif not _all_downstream_outputs_exist(dep):
                return False
        return True

    @lru_cache(maxsize=None)
    def _should_run(target):
        if spec_hashes.has_changed(target) is not None:
            logger.debug("Target %s has a changed spec", target)
            return True

        # If I have no outputs, but I have inputs, I should probably only run if my input
        # changed, but I don't have any output files to compare with, so I'll just run
        # every time.
        if not target.outputs:
            logger.debug(
                "Target %s has no outputs and will always be scheduled", target
            )
            return True

        all_outputs = set(target.flattened_outputs())
        non_temp_outputs = all_outputs - graph.temporary
        temp_outputs = all_outputs & graph.temporary

        # Targets with non-temporary outputs
        if non_temp_outputs:
            missing_output = next(
                (p for p in non_temp_outputs if not fs.exists(p)), None
            )
            if missing_output is not None:
                logger.debug(
                    "Target %s is missing output file %s", target, missing_output
                )
                return True

            non_temp_inputs = set(target.flattened_inputs()) - graph.temporary
            missing_input = next((p for p in non_temp_inputs if not fs.exists(p)), None)
            if missing_input is not None:
                logger.debug(
                    "Target %s is missing input file %s", target, missing_input
                )
                return True

            if _is_outdated(non_temp_inputs, non_temp_outputs):
                logger.debug("Target %s is not up-to-date", target)
                return True

            logger.debug("Target %s is up-to-date", target)
            return False

        # Targets with only temporary outputs
        missing_temp = [p for p in temp_outputs if not fs.exists(p)]
        if missing_temp:
            if _all_downstream_outputs_exist(target):
                logger.debug(
                    "Target %s has missing temporary outputs but all downstream outputs exist",
                    target,
                )
                return False

            logger.debug(
                "Target %s is missing temporary output file %s",
                target,
                missing_temp[0],
            )
            return True

        existing_inputs = [
            p for p in target.flattened_inputs() if p in graph.temporary or fs.exists(p)
        ]
        if _is_outdated(existing_inputs, temp_outputs):
            logger.debug("Target %s (temp outputs) is not up-to-date", target)
            return True

        logger.debug("Target %s is up-to-date", target)
        return False

    def _schedule(target):
        submitted_deps = []

        if not no_deps:
            for dep in sorted(graph.dependencies[target], key=lambda t: t.name):
                status = _cached_schedule(dep)
                if status in SUBMITTED_STATES:
                    submitted_deps.append(dep)

        if force:
            logger.debug("Target %s is being forcibly submitted", target)
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        if status_func(target) == BackendStatus.SUBMITTED:
            logger.debug("Target %s is already submitted", target)
            return Status.SUBMITTED

        if status_func(target) == BackendStatus.RUNNING:
            logger.debug("Target %s is already running", target)
            return Status.RUNNING

        if status_func(target) == BackendStatus.FAILED:
            submit_func(target, dependencies=submitted_deps)
            return Status.FAILED

        if status_func(target) == BackendStatus.CANCELLED:
            submit_func(target, dependencies=submitted_deps)
            return Status.CANCELLED

        if submitted_deps:
            logger.debug(
                "Target %s will be submitted because of dependency %s",
                target,
                submitted_deps[0],
            )
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        if _should_run(target):
            submit_func(target, dependencies=submitted_deps)
            return Status.SHOULDRUN

        return Status.COMPLETED

    def _cached_schedule(target):
        if target not in schedule_cache:
            schedule_cache[target] = _schedule(target)
        return schedule_cache[target]

    for target in sorted(endpoints, key=lambda t: t.name):
        _cached_schedule(target)

    return schedule_cache


def _submit_dryrun(target, dependencies, backend, spec_hashes):
    logger.info("Would submit %s", target)


def _submit_noop(target, dependencies, backend, spec_hashes):
    pass


def submit_backend(target, dependencies, backend, spec_hashes):
    """Prepare and submit `target` with `dependencies` to `backend`.

    Will prepare the target for submission by injecting option defaults from the
    backend, check for unsupported options, and removing options with a `None`
    value.

    This is the primary way to submit a target. Do not call :func:`submit`
    directly on the backend, unless you want to manually deal with with
    injection of option defaults.
    """

    # logger.info("Submitting target %s", target)

    new_options = {}
    if hasattr(backend, "target_defaults"):
        new_options = dict(backend.target_defaults)
    new_options.update(target.options)

    for option_name, option_value in list(new_options.items()):
        if option_name not in backend.target_defaults.keys():
            logger.warning(
                "Option '%s' used in '%s' is not supported by backend. Ignored.",
                option_name,
                target.name,
            )
            del new_options[option_name]
        elif option_value is None:
            del new_options[option_name]
    target.options = new_options

    backend.submit(target, dependencies)
    spec_hashes.update(target)

    if hasattr(backend, "get_tracked_id"):
        logger.info(
            "Submitted target %s (id: %s)",
            target,
            backend.get_tracked_id(target),
        )
    else:
        logger.info("Submitted target %s", target)


def submit_workflow(
    endpoints,
    graph,
    fs,
    spec_hashes,
    backend,
    dry_run=False,
    force=False,
    no_deps=False,
):
    """Submit a workflow to a backend."""
    submit_func = partial(
        _submit_dryrun if dry_run else submit_backend,
        backend=backend,
        spec_hashes=spec_hashes,
    )
    schedule(
        endpoints,
        graph,
        fs,
        spec_hashes,
        status_func=backend.status,
        submit_func=submit_func,
        force=force,
        no_deps=no_deps,
    )


def get_status_map(graph, fs, spec_hashes, backend, endpoints=None):
    """Get the status of each targets in the graph."""
    submit_func = partial(_submit_noop, backend=backend, spec_hashes=spec_hashes)
    return schedule(
        endpoints if endpoints is not None else graph.endpoints(),
        graph,
        fs,
        spec_hashes,
        status_func=backend.status,
        submit_func=submit_func,
    )
