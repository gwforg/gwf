#!/usr/bin/env bash
#
# Upload conda packages for the current Python version and all supported platforms (osx, win, linux).
#
# Author:
#   Dan Søndergaard <das@birc.au.dk>
#

set -eu
set -o pipefail

if [[ $TRAVIS_BRANCH = "master" ]]; then
    LABEL="main"
else
    LABEL="dev"
fi

conda convert --platform all $HOME/miniconda/conda-bld/*/*.tar.bz2 -o $HOME/miniconda/conda-bld/
anaconda -t $ANACONDA_TOKEN upload --label $LABEL --user gwforg $HOME/miniconda/conda-bld/*/*.tar.bz2
