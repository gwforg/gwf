#!/usr/bin/env bash

VERSION=`python setup.py --version`

# Setup Miniconda...
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r

# Configure and update...
conda config --set always_yes yes --set changeps1 no
conda config --add channels dansondergaard
conda update -q conda
conda install anaconda-client

cd /tmp

# Build skeleton from PyPI package...
conda skeleton pypi --pypi-url https://testpypi.python.org/pypi gwf --version $VERSION
conda build --python 3.4 gwf/
conda convert --platform all $HOME/miniconda/conda-bld/linux-64/gwf-$VERSION-py34_0.tar.bz2

for f in `ls */*.tar.bz2`; do
    anaconda upload $f
done