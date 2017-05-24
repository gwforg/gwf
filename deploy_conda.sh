#!/usr/bin/env bash
#
# Upload conda packages for the current Python version and all supported platforms (osx, win, linux).
#
# Author:
#   Dan Søndergaard <das@birc.au.dk>
#

set -eu
set -o pipefail

conda convert --platform all $HOME/conda-bld/linux-64/gwf-*-*_0.tar.bz2 -o conda-bld
anaconda -u gwforg -t $ANACONDA_TOKEN upload conda-bld/*/*.tar.bz2
