#!/bin/bash

# set desired vector length
N=20000000

# configure build
cmake .

# build executables
make

# run serial test
echo "  "
echo "running serial version, dot_prod.serial:"
./dot_prod.serial $N

# run CUDA test
echo "  "
echo "running CUDA version, dot_prod.cuda:"
./dot_prod.cuda $N
