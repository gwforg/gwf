import nox

import os.path
import shutil


def setup_test_env(session):
    session.install("--upgrade", "pip", "setuptools")
    session.install("-r", "requirements.txt", "-r", "requirements-dev.txt")
    session.run("pip", "install", "--no-deps", "--editable", ".")


@nox.session(python=["3.5", "3.6", "3.7", "3.8"])
def test(session):
    setup_test_env(session)

    session.run(
        "pytest",
        "--doctest-modules",
        "--cov-config=.coveragerc",
        "--cov=src/gwf",
        "tests/",
    )


@nox.session(python="3.6")
def black(session):
    session.install("black")
    session.run("black", "--check", "--diff", "src/", "tests/")


@nox.session
def lint(session):
    setup_test_env(session)
    session.run("flake8", "src/gwf")


@nox.session(python="3.6")
def docs(session):
    setup_test_env(session)
    session.chdir("docs")
    if os.path.exists("_build"):
        shutil.rmtree("_build")
    session.run("sphinx-build", "-b", "dirhtml", ".", "_build/dirhtml")
