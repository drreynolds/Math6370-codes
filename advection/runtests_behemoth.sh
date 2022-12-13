#!/bin/bash

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.7.00/gcc-11.3.0/lib/cmake/Kokkos

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

# run OpenMP test with 16 threads
echo "  "
echo "running OpenMP version, advection.openmp, with 16 threads:"
OMP_NUM_THREADS=16 ./advection.openmp
