#!/bin/bash

# set desired number of intervals
N=10000000

# set desired numbers of OpenMP threads
THREADS=(1 2 4 8 16 32)

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.6.00/gcc-11.2.0/lib/cmake/Kokkos

# build executables
make

# run serial test
echo "  "
echo "running serial version, pi_comp.serial:"
./pi_comp.serial $N

# run OpenMP tests
for t in "${THREADS[@]}"
do
    echo "  "
    echo "running OpenMP version, pi_comp.openmp, with $t threads:"
    OMP_NUM_THREADS=$t ./pi_comp.openmp $N

done

# run CUDA test
echo "  "
echo "running CUDA version, pi_comp.cuda:"
./pi_comp.cuda $N
