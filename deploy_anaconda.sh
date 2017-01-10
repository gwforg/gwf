#!/bin/bash

set -e

echo "Converting conda package..."
conda convert --platform all $HOME/miniconda3/conda-bld/linux-64/gwf-*.tar.bz2 --output-dir conda-bld/

echo "Deploying to Anaconda.org..."
anaconda -t $ANACONDA_TOKEN upload conda-bld/**/PACKAGENAME-*.tar.bz2

echo "Successfully deployed to Anaconda.org."
exit 0