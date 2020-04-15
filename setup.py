from setuptools import setup, find_packages

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name="gwf",
    version="1.7.2",
    author="Dan Søndergaard",
    author_email="das@birc.au.dk",
    description="A flexible, pragmatic workflow tool.",
    long_description=long_description,
    url="https://gwf.app/",
    project_urls={
        "Issues": "https://github.com/gwforg/gwf/issues",
        "Source code": "https://github.com/gwforg/gwf",
    },
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.5",
    install_requires=["click", "click-plugins"],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Operating System :: OS Independent",
        "Topic :: Utilities",
        "Topic :: System :: Distributed Computing",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
    entry_points={
        "console_scripts": ["gwf = gwf.cli:main"],
        "gwf.plugins": [
            "init = gwf.plugins.init:init",
            "run = gwf.plugins.run:run",
            "config = gwf.plugins.config:config",
            "status = gwf.plugins.status:status",
            "info = gwf.plugins.info:info",
            "logs = gwf.plugins.logs:logs",
            "clean = gwf.plugins.clean:clean",
            "workers = gwf.plugins.workers:workers",
            "cancel = gwf.plugins.cancel:cancel",
            "touch = gwf.plugins.touch:touch",
        ],
        "gwf.backends": [
            "slurm = gwf.backends.slurm:SlurmBackend",
            "sge = gwf.backends.sge:SGEBackend",
            "local = gwf.backends.local:LocalBackend",
            "testing = gwf.backends.testing:TestingBackend",
        ],
    },
)
