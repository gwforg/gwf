import nox

import os.path
import shutil


@nox.session(python=["3.8", "3.9", "3.10"])
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

    session.notify("coverage")


# @nox.session(python="3.10")
# def coverage(session):
#     session.install("coverage[toml]")
#     session.run("coverage", "combine")
#     session.run("coverage", "report", "--fail-under=65", "--show-missing")
#     session.run("coverage", "erase")

#     token = os.getenv("COVERALLS_REPO_TOKEN")
#     if token:
#         session.install("coveralls")
#         session.run("coveralls", env={"COVERALLS_REPO_TOKEN": token})


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
    session.run("flake8", *files)


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
    sphinx_args = ["-b", "html", ".", "_build/html"]

    if not session.interactive:
        sphinx_cmd = "sphinx-build"
    else:
        sphinx_cmd = "sphinx-autobuild"
        sphinx_args.insert(0, "--open-browser")

    session.run(sphinx_cmd, *sphinx_args)
