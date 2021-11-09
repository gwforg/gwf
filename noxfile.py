import nox

import os.path
import shutil


@nox.session(python=["3.8", "3.9", "3.10"])
def test(session):
    session.install("--upgrade", "pip", "flit")
    session.run("flit", "install", "--deps", "production")

    session.run(
        "pytest",
        "--cov-config",
        "pyproject.toml",
        "--cov",
        "src/gwf",
        "tests/",
        env={"COVERAGE_FILE": f".coverage.{session.python}"},
    )

    session.notify("coverage")


@nox.session(python="3.10")
def coverage(session):
    session.install("--upgrade", "pip", "flit")
    session.run("flit", "install", "--deps", "production", "--extras", "test")
    session.run("coverage", "combine")
    session.run("coverage", "report", "--fail-under=95", "--show-missing")
    session.run("coverage", "erase")


@nox.session
def format(session):
    session.install("--upgrade", "pip", "flit")
    session.run("flit", "install", "--deps", "production", "--extras", "dev")

    files = ["src/gwf", "tests"]
    session.run("black", *files)
    session.run("isort", *files)


@nox.session
def lint(session):
    session.install("--upgrade", "pip", "flit")
    session.run("flit", "install", "--deps", "production", "--extras", "dev,test")

    files = ["src/gwf", "tests"]
    session.run("black", "--check", *files)
    session.run("isort", "--check", *files)
    session.run("flake8", *files)


@nox.session(python="3.10")
def docs(session):
    session.install("--upgrade", "pip", "flit")
    session.run("flit", "install", "--deps", "production", "--extras", "docs")

    session.cd("docs")
    if os.path.exists("_build"):
        shutil.rmtree("_build")

    sphinx_args = ["-b", "html", "-W", ".", "_build/html"]

    if not session.interactive:
        sphinx_cmd = "sphinx-build"
    else:
        sphinx_cmd = "sphinx-autobuild"
        sphinx_args.insert(0, "--open-browser")

    session.run(sphinx_cmd, *sphinx_args)
