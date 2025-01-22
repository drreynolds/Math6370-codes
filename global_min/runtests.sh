#!/bin/bash

# configure build
cmake .

# build executables
make

# run serial test
echo "  "
echo "running serial version, glob_min.serial:"
./glob_min.serial

# run CUDA test
echo "  "
echo "running CUDA version, glob_min.cuda:"
./glob_min.cuda
