import ez_setup
ez_setup.use_setuptools()

import os
from setuptools import setup, find_packages

# Utility function to read the README file.  Used for the
# long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put
# a raw string in below ...

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gwf",
    version = "0.7.0",

    packages = find_packages(where='src'),
    package_dir = {'': 'src'},
    scripts=['scripts/gwf', 'scripts/gwf-config'],

    #test_suite='tests',

    # metadata for upload to PyPI
    author = "Thomas Mailund",
    author_email = "mailund@birc.au.dk",
    license = "GPLv3",
    keywords = "grid computing workflow",
    url = "https://mailund.github.io/gwf",
    description = ("Grid WorkFlow - a make-like system for computer grids."),
    long_description = read('README'),

    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Topic :: Utilities",
        "Topic :: System :: Distributed Computing",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPLv3)",
        "Programming Language :: Python",
    ],
)
