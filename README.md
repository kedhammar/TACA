<p align="center">
  <a href="https://github.com/SciLifeLab/TACA">
    <img width="512" height="175" src="artwork/logo.png"/>
  </a>
</p>

## Tool for the Automation of Cleanup and Analyses

[![PyPI version](https://badge.fury.io/py/taca.svg)](http://badge.fury.io/py/taca)
[![Build Status](https://travis-ci.org/SciLifeLab/TACA.svg?branch=master)](https://travis-ci.org/SciLifeLab/TACA)
[![Documentation Status](https://readthedocs.org/projects/taca/badge/?version=latest)](https://readthedocs.org/projects/taca/?badge=latest)
[![codecov](https://codecov.io/gh/scilifelab/taca/branch/master/graph/badge.svg)](https://codecov.io/gh/scilifelab/taca)

This package contains several tools for projects and data management in the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) in Stockholm, Sweden.

### Install for development
You can install your own fork of taca in for instance a local conda environment for development. Provided you have conda installed:

```
# clone the repo
git clone https://github.com/<username>/TACA.git

# create an environment
conda create -n taca_dev python=2.7
conda activate taca_dev

# install TACA and dependencies for developoment
cd TACA
python setup.py develop
pip install -r ./requirements-dev.txt 

# Check that tests pass:
cd tests && nosetests -v -s
```

There is also a [plugin for the deliver command](https://github.com/SciLifeLab/taca-ngi-pipeline). To install this in the same development environment:

```
# Install taca delivery plugin for development
git clone https://github.com/<username>/TACA.git
cd ../taca-ngi-pipeline
python setup.py develop
pip install -r ./requirements-dev.txt

# add required config files and env for taca delivery plugin
echo "foo:bar" >> ~/.ngipipeline/ngi_config.yaml 
mkdir ~/.taca && cp tests/data/taca_test_cfg.yaml ~/.taca/taca.yaml
export CHARON_BASE_URL="http://tracking.database.org"
export CHARON_API_TOKEN="charonapitokengoeshere"

# Check that tests pass:
cd tests && nosetests -v -s
```

For a more detailed documentation please go to [the documentation page](http://taca.readthedocs.org/en/latest/).
