#!/bin/bash

# set desired vector length
N=10000000

# configure build
cmake .

# build executables
make

# run serial test
echo "  "
echo "running serial version, axpy.serial:"
./axpy.serial $N

# run CUDA test
echo "  "
echo "running CUDA version, axpy.cuda:"
./axpy.cuda $N
