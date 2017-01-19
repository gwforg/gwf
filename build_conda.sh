#!/usr/bin/env bash
#
# Builds and uploads conda packages for the current Python version and
# all supported platforms (osx, win, linux).

set -e

cd /tmp

# Setup Miniconda...
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /tmp/miniconda
export PATH="/tmp/miniconda/bin:$PATH"

# Configure and update...
conda config --set always_yes yes --set changeps1 no
conda config --add channels dansondergaard
conda update -q conda
conda install anaconda-client conda-build

# Build skeleton from PyPI package...
conda skeleton pypi --pypi-url $PYPI_URL gwf
conda build --python $TRAVIS_PYTHON_VERSION gwf/
mkdir /tmp/builds
conda convert --platform all /tmp/miniconda/conda-bld/linux-64/gwf-*-*_0.tar.bz2 -o /tmp/builds

for f in `ls /tmp/builds/*/*.tar.bz2`; do
    anaconda -t $ANACONDA_TOKEN upload $f
done