#!/usr/bin/env bash
#
# Builds and uploads conda packages for the current Python version and
# all supported platforms (osx, win, linux).

set -e

conda convert --platform all $HOME/conda-bld/linux-64/gwf-*-*_0.tar.bz2 -o conda-bld
anaconda -u gwforg -t $ANACONDA_TOKEN upload conda-bld/*/*.tar.bz2
