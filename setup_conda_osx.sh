#!/bin/bash

# Workaround for https://github.com/travis-ci/travis-ci/issues/6307, which
# caused the following error on MacOS X workers:
#
# /Users/travis/build.sh: line 109: shell_session_update: command not found
#
rvm get head

# Install conda
# http://conda.pydata.org/docs/travis.html#the-travis-yml-file

if ! type conda > /dev/null; then
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $HOME/cache/miniconda
    ls -s $HOME/cache/miniconda $HOME/miniconda
fi
export PATH="$HOME/miniconda/bin:$PATH"

# Install common Python dependencies
source "$( dirname "${BASH_SOURCE[0]}" )"/setup_dependencies_common.sh