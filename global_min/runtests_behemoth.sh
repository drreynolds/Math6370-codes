#!/bin/bash

# set desired numbers of OpenMP threads
THREADS=(1 2 4 8 16 32)

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.6.00/gcc-11.2.0/lib/cmake/Kokkos

# build executables
make

# run serial test
echo "  "
echo "running serial version, glob_min.serial:"
./glob_min.serial

# run OpenMP tests
for t in "${THREADS[@]}"
do
    echo "  "
    echo "running OpenMP version, glob_min.openmp, with $t threads:"
    OMP_NUM_THREADS=$t ./glob_min.openmp

done

# run CUDA test
echo "  "
echo "running CUDA version, glob_min.cuda:"
./glob_min.cuda
