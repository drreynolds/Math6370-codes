#!/bin/bash

# set desired number of intervals
N=10000000

# configure build
cmake .

# build executables
make

# run serial test
echo "  "
echo "running serial version, pi_comp.serial:"
./pi_comp.serial $N

# run CUDA test
echo "  "
echo "running CUDA version, pi_comp.cuda:"
./pi_comp.cuda $N
