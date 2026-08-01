import asyncio
import logging
import os
import signal
import stat
import sys
import tempfile
import time
from gwf.executors import deserialize
from gwf.exec.isolation import IsolationError, execution_context, isolate_network


logger = logging.getLogger(__name__)


async def forward(src, dst, bufsize=2**16, flush_sec=10):
    """Forward `src` to `dst` with buffering and time-based flushing.

    This will forward `src` to `dst`, but buffer writes to `dst` in-memory. Data
    will be flushed to `dst` when the buffer is filled or after `flush_sec`
    seconds.
    """
    buf = bytearray(bufsize)
    buf_used = 0

    last_flush = time.monotonic()
    done = False
    while not done:
        try:
            data = await asyncio.wait_for(src.read(bufsize - buf_used), 1)
            if not data:
                done = True
            buf[buf_used : buf_used + len(data)] = data
            buf_used += len(data)
        except asyncio.TimeoutError:
            pass

        now = time.monotonic()
        if done or buf_used == bufsize or now - last_flush > flush_sec:
            if buf_used:
                dst.write(buf[:buf_used])
                dst.flush()
            last_flush = now
            buf_used = 0


async def execute_command(cmd, working_dir, environ=None):
    proc = await asyncio.create_subprocess_exec(
        *cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
        bufsize=0,
        env=environ,
        cwd=working_dir,
        start_new_session=True,
    )

    stdout_task = asyncio.create_task(forward(proc.stdout, sys.stdout.buffer))
    stderr_task = asyncio.create_task(forward(proc.stderr, sys.stderr.buffer))
    proc_task = proc.wait()

    try:
        await asyncio.gather(proc_task, stdout_task, stderr_task)
    except BaseException:
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        if proc.returncode is None:
            try:
                await asyncio.wait_for(proc.wait(), timeout=5)
            except asyncio.TimeoutError:
                pass
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        if proc.returncode is None:
            await proc.wait()
        await asyncio.gather(stdout_task, stderr_task, return_exceptions=True)
        raise
    return proc.returncode


def execute_script(script_file, workflow_root, debug_mode):
    target = deserialize(script_file)
    logger.info("executing target %r", target)
    return_code = 0
    with execution_context(target) as execution:
        spec_file, spec_path = tempfile.mkstemp(
            prefix="gwf_",
            dir=execution.spec_directory,
            text=True,
        )

        logger.debug("writing script with spec to %s", spec_path)
        try:
            spec_file = os.fdopen(spec_file, "w")
            spec_file.write("#!/bin/bash\n\n")
            spec_file.write("set -e\n")

            # TODO: At some point introduce a 'strict mode' setting which enables these.
            # spec_file.write("set -u\n")
            # spec_file.write("set -o pipefail\n")

            if debug_mode:
                spec_file.write("set -x\n")

            spec_file.write("\n")
            spec_file.write(target.spec)
            spec_file.flush()
            spec_file.close()
            os.chmod(spec_path, stat.S_IRUSR | stat.S_IEXEC)

            if debug_mode:
                logger.debug(
                    "will execute the following script:\n%s",
                    open(spec_path).read(),
                )

            if execution.isolated:
                cmd = target.executor.get_isolated_command(
                    spec_path,
                    workflow_root,
                    target.working_dir,
                )
            else:
                cmd = target.executor.get_command(spec_path, workflow_root)

            if execution.isolated and not target.isolation.network:
                cmd = isolate_network(cmd, target.name)

            environ = os.environ.copy()
            environ["GWF_TARGET_NAME"] = target.name

            logger.debug("running command %s", cmd)
            # Filesystem restrictions belong at this child process boundary.
            # The gwf.exec parent must retain access to stage inputs and publish
            # outputs.
            return_code = asyncio.run(
                execute_command(cmd, execution.working_dir, environ)
            )
            logger.debug("command finished with return code %s", return_code)
            if return_code == 0:
                execution.publish_outputs()
        finally:
            logger.debug("removing spec file at %s", spec_path)
            os.remove(spec_path)

    return return_code


def main():
    debug_mode = os.getenv("GWF_EXEC_DEBUG_MODE", False)
    workflow_root = os.environ["GWF_EXEC_WORKFLOW_ROOT"]

    logging.basicConfig(level=logging.DEBUG if debug_mode else logging.WARNING)
    logger.info("this is gwf-exec, hello!")

    def terminate(signum, frame):
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, terminate)

    script_path = sys.argv[-1]
    with open(script_path) as script_file:
        try:
            return_code = execute_script(script_file, workflow_root, debug_mode)
        except IsolationError as error:
            if debug_mode:
                raise
            error.show()
            sys.exit(error.exit_code)
    sys.exit(return_code)


main()
