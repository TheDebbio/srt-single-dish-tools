#!/bin/bash

# Install conda
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file

# This if accounts for a cached installation of anaconda

if ! type conda > /dev/null; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/cache/miniconda
    ls -s $HOME/cache/miniconda $HOME/miniconda
fi
export PATH="$HOME/miniconda/bin:$PATH"

# Install common Python dependencies
source "$( dirname "${BASH_SOURCE[0]}" )"/setup_dependencies_common.sh

if [[ $SETUP_XVFB == True ]]; then
    export DISPLAY=:99.0
    sh -e /etc/init.d/xvfb start
fi
