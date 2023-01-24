import nox

import os
import shutil
from fnmatch import fnmatch


PURGE_PATTERNS = [
    "*/__pycache__/*",
    "docs/_build/*",
    ".coverage/*",
    "*/.egg-info/*",
    "*/.egg/*",
    "*/.eggs/*",
    "*/.gwf/*",
    "*/.gwfconf.json",
    "*/.pytest_cache/*",
    "build/*",
    "conda-bld/*",
    "dist/*",
    "*.pyc",
]


def matches_pattern(path, patterns):
    for pattern in patterns:
        if fnmatch(path, pattern):
            return True
    return False


@nox.session
def clean(session):
    def _delete(path, patterns):
        for entry in os.scandir(path):
            if entry.is_dir():
                _delete(entry.path, patterns)
            elif entry.is_file() and matches_pattern(entry.path, patterns):
                os.remove(entry.path)
        if matches_pattern(path, patterns):
            os.rmdir(path)

    _delete(".", PURGE_PATTERNS)


@nox.session
def build(session):
    session.install("flit")
    session.run("flit", "build", "--no-setup-py", "--format", "wheel")


@nox.session
def package(session):
    session.run("conda", "build", "--output-folder", "conda-bld/", "conda/")


@nox.session
def publish(session):
    session.run(
        "anaconda",
        "-t",
        "${ANACONDA_TOKEN}",
        "upload",
        "--force",
        "--no-progress",
        "--user",
        "gwforg",
        "conda-bld/*/*.tar.bz2",
    )


@nox.session(python=["3.7", "3.8", "3.9", "3.10"])
def test(session):
    session.install("flit")
    session.run("flit", "install", "-s", "--deps", "production")

    session.install(
        "flake8",
        "pytest",
        "pytest-click",
        "pytest-cov",
        "pytest-flake8",
        "pytest-mock",
        "pytest-runner",
    )
    session.warn(session.virtualenv)
    session.run(
        "pytest",
        "--cov-config",
        "pyproject.toml",
        "--cov",
        "gwf",
        "tests/",
        env={"COVERAGE_FILE": f".coverage.{session.python}"},
    )

    if session.interactive:
        session.notify("coverage")


@nox.session(python="3.10")
def coverage(session):
    session.install("coverage[toml]")
    session.run("coverage", "combine")
    session.run("coverage", "report", "--fail-under=65", "--show-missing")
    session.run("coverage", "erase")

    token = os.getenv("COVERALLS_REPO_TOKEN")
    if token:
        session.install("coveralls")
        session.run("coveralls", env={"COVERALLS_REPO_TOKEN": token})


@nox.session
def format(session):
    session.install("black", "isort")

    files = ["src/gwf", "tests"]
    session.run("black", *files)
    session.run("isort", *files)


@nox.session
def lint(session):
    session.install("black", "isort", "flake8")

    files = ["src/gwf", "tests"]
    session.run("black", "--check", *files)
    session.run("isort", "--check", *files)

    # TODO: Also lint tests at some point...
    session.run("flake8", "src/gwf")


@nox.session(python="3.10")
def docs(session):
    session.install("flit")
    session.run("flit", "install", "-s", "--deps", "production")

    session.install(
        "sphinx",
        "sphinx-autobuild",
        "furo",
    )

    session.cd("docs")
    if os.path.exists("_build"):
        shutil.rmtree("_build")

    # TODO: Add -W to treat warnings as errors
    sphinx_args = ["-b", "dirhtml", ".", "_build/dirhtml"]

    if not session.interactive:
        sphinx_cmd = "sphinx-build"
    else:
        sphinx_cmd = "sphinx-autobuild"
        sphinx_args.insert(0, "--open-browser")

    session.run(sphinx_cmd, *sphinx_args)
