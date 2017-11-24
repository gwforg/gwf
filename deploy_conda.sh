#!/usr/bin/env bash
#
# Upload conda packages for the current Python version and all supported platforms (osx, win, linux).
#
# Author:
#   Dan Søndergaard <das@birc.au.dk>
#

set -eu
set -o pipefail

conda convert --platform all $HOME/miniconda/conda-bld/*/*.tar.bz2 -o $HOME/miniconda/conda-bld/
anaconda -t $ANACONDA_TOKEN upload --user gwforg $HOME/miniconda/conda-bld/*/*.tar.bz2
