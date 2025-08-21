import hashlib
from pathlib import Path
from typing import Tuple

import click


def get_logs_dir(working_dir: Path):
    return working_dir.joinpath(".gwf", "logs")


def _get_log_path_without_extension(working_dir: Path, name: str):
    digest = hashlib.sha256(name.encode("utf-8")).hexdigest()
    suffix = f"{digest[0]}/{digest[1]}/{name}"
    return get_logs_dir(working_dir).joinpath(suffix)


def get_log_paths(working_dir: Path | str, name: str) -> Tuple[Path, Path]:
    base_path = _get_log_path_without_extension(Path(working_dir), name)
    stdout_path = base_path.with_suffix(".stdout")
    stderr_path = base_path.with_suffix(".stderr")
    return stdout_path, stderr_path


def prepare_log_storage_for_target(working_dir: Path, name: str):
    path, _ = get_log_paths(working_dir, name)
    path.parent.mkdir(parents=True, exist_ok=True)


def init_log_storage(working_dir: Path | str):
    working_dir = Path(working_dir)
    logs_dir = get_logs_dir(working_dir)
    logs_dir.mkdir(exist_ok=True, parents=True)

    # Migrate log files from the old, flat hierarchy to the new, hash-based
    # hierarcy.
    migrate_files = []
    for entry in logs_dir.iterdir():
        if entry.suffix in (".stdout", ".stderr"):
            migrate_files.append(entry)

    if migrate_files:
        with click.progressbar(migrate_files, label="Migrating log files") as bar:
            for entry in bar:
                new_path = prepare_log_storage_for_target(working_dir, entry.stem)
                new_path = _get_log_path_without_extension(
                    working_dir, entry.stem
                ).with_suffix(entry.suffix)
                prepare_log_storage_for_target(working_dir, entry.stem)
                entry.rename(new_path)


def clean_logs(working_dir: Path | str, graph):
    working_dir = Path(working_dir)

    active_paths = set()
    for name in graph.targets.keys():
        out, err = get_log_paths(working_dir, name)
        active_paths.add(out)
        active_paths.add(err)

    emptied_dirs = set()
    for root, dirs, files in get_logs_dir(working_dir).walk(top_down=False):
        files_removed = 0
        dirs_removed = 0

        for name in dirs:
            path = root.joinpath(name)
            if path in emptied_dirs:
                path.rmdir()
                dirs_removed += 1

        for name in files:
            path = root.joinpath(name)
            if path not in active_paths:
                path.unlink()
                files_removed += 1

        if files_removed == len(files) and dirs_removed == len(dirs):
            emptied_dirs.add(root)
