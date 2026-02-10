import logging
from functools import partial

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
DEPENDENT_SHOULDRUN_STATES = {BackendStatus.FAILED, BackendStatus.CANCELLED}


def should_run(target, graph, fs, spec_hashes):
    if spec_hashes.has_changed(target) is not None:
        logger.debug("Target %s has a changed spec", target)
        return True

    # If I have no outputs, but I have inputs, I should probably only run if my input
    # changed, but I don't have any output files to compare with, so I'll just run
    # every time.
    if not target.outputs:
        logger.debug("Target %s has no outputs and will always be scheduled", target)
        return True

    # Only check non-temp outputs for existence
    non_temp_outputs = set(target.flattened_outputs()) - graph.temporary
    for path in non_temp_outputs:
        if not fs.exists(path):
            logger.debug("Target %s is missing output file %s", target, path)
            return True

    # Only consider non-temp inputs for timestamp comparison
    non_temp_inputs = set(target.flattened_inputs()) - graph.temporary
    youngest_in_ts, _ = max(
        ((fs.changed_at(path), path) for path in non_temp_inputs),
        default=(float("-inf"), None),
    )

    oldest_out_ts, _ = min(
        ((fs.changed_at(path), path) for path in non_temp_outputs),
        default=(float("inf"), None),
    )

    if youngest_in_ts > oldest_out_ts:
        logger.debug("Target %s is not up-to-date", target)
        return True

    logger.debug("Target %s is up-to-date", target)
    return False


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
    should_run_cache = {}
    schedule_cache = {}

    def _cached_should_run(target):
        if target not in should_run_cache:
            should_run_cache[target] = should_run(target, graph, fs, spec_hashes)
        return should_run_cache[target]

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

        # For targets that produce temp files, check if any dependent needs to run
        if target.temp:
            for dependent in graph.dependents[target]:
                if status_func(
                    dependent
                ) in DEPENDENT_SHOULDRUN_STATES or _cached_should_run(dependent):
                    logger.debug(
                        "Target %s will be submitted because dependent %s needs it",
                        target,
                        dependent,
                    )
                    submit_func(target, dependencies=submitted_deps)
                    return Status.SHOULDRUN

        if _cached_should_run(target):
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
