#!/bin/bash

# configure build
cmake .

# build executables
make

# run CUDA test
echo "  "
echo "running CUDA version, advection.cuda:"
./advection.cuda

# run hybrid CUDA+OpenMP test
echo "  "
echo "running hybrid CUDA+OpenMP version, advection.hybrid:"
./advection.hybrid
