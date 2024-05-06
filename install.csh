#!/usr/bin/env bash

module load mpi-x84
conda create -n amuse
conda activate amuse
mkdir amuse
cd amuse
git clone https://github.com/amusecode/amuse.git .
conda install python
conda install numpy
conda install scipy
conda install h5py
conda install pytest
conda install docutils
conda install mpfr

./configure --prefix=$CONDA_PREFIX
make DOWNLOAD_CODES=1

pip install -e .
python setup.py develop_build
