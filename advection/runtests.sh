#!/bin/bash

# configure build
cmake .

# build executables
make

# run CUDA test
echo "  "
echo "running CUDA version, advection.cuda:"
./advection.cuda
