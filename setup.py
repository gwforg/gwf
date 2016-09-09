# -*- coding: utf-8 -*-

from __future__ import (absolute_import, print_function, division,
                        unicode_literals)

import sys
import os

from setuptools import find_packages, setup


requires = ['six>=1.10']
if sys.version_info.major == 2:
    requires.append('contextdecorator')

# Utility function to read the README file.  Used for the
# long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put
# a raw string in below ...


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="gwf",
    version="1.0.0",

    packages=find_packages(),
    scripts=['scripts/gwf'],

    test_suite='tests',

    # metadata for upload to PyPI
    author="Thomas Mailund, Dan Søndergaard",
    author_email="mailund@birc.au.dk, das@birc.au.dk",
    license="GPLv3",
    keywords="grid computing workflow",
    url="https://mailund.github.io/gwf",
    description="A flexible, pragmatic workflow tool.",
    long_description=read('README.rst'),
    install_requires=requires,
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
