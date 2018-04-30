#!/usr/bin/env bash
#
# Upload conda packages for the current Python version and all supported platforms (osx, win, linux).
#
# Author:
#   Dan Søndergaard <das@birc.au.dk>
#

set -eu
set -o pipefail

wget -q https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
chmod +x miniconda.sh
./miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"

conda update --yes -q conda
conda config --set always_yes true
conda config --set anaconda_upload no
conda config --add channels gwforg

conda install -q python=$TRAVIS_PYTHON_VERSION pip=10.0.1 conda-build=3.0.15 anaconda-client=1.6.0
conda build --python $TRAVIS_PYTHON_VERSION conda/
conda convert --platform all $HOME/miniconda/conda-bld/*/*.tar.bz2 -o $HOME/miniconda/conda-bld/
anaconda -t $ANACONDA_TOKEN upload --user gwforg $HOME/miniconda/conda-bld/*/*.tar.bz2
