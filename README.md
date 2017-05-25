# SRT Single dish tools #
![status](https://gitlab.com/matteobachetti/srt-single-dish-tools/badges/master/build.svg)
![coverage](https://gitlab.com/matteobachetti/srt-single-dish-tools/badges/master/coverage.svg)

## Installation

### Anaconda and virtual environment (recommended but optional)

We strongly suggest to install the
[Anaconda](https://www.continuum.io/downloads) Python distribution.
Once the installation has finished, you should have a working `conda`
command in your shell. First of all, create a new environment:

    $ conda create -n py3 python=3

load the new environment:

    $ source activate py3

and install the dependencies:

    (py3) $ conda install matplotlib h5py astropy scipy numpy

### Cloning and installation

Clone the repository:

    (py3) $ cd /my/software/directory/
    (py3) $ git clone https://matteobachetti@gitlab.com/matteobachetti/srt-single-dish-tools.git

or if you have deployed your SSH key to Gitlab:

    (py3) $ git git@gitlab.com:matteobachetti/srt-single-dish-tools.git

Then:

    (py3) $ cd srt-single-dish-tools
    (py3) $ python setup.py install

That's it. After installation has ended, you can verify that software is
installed by executing:

    (py3) $ SDTimage -h

If the help message appears, you're done!

### Updating

To update the code, simply run `git pull` and reinstall:

    (py3) $ git pull
    (py3) $ python setup.py install

### Contribution guidelines ###

[Why writing contribution guidelines, when Astropy has made such a good job already?](http://docs.astropy.org/en/stable/development/codeguide.html)
