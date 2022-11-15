#!/bin/bash

# set desired vector length
N=10000000

# set desired numbers of OpenMP threads
THREADS=(1 2 4 8 16 32)

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.6.00/gcc-11.2.0/lib/cmake/Kokkos

# build executables
make

# run serial test
echo "  "
echo "running serial version, axpy.serial:"
./axpy.serial $N

# run OpenMP tests
for t in "${THREADS[@]}"
do
    echo "  "
    echo "running OpenMP version, axpy.openmp, with $t threads:"
    OMP_NUM_THREADS=$t ./axpy.openmp $N

done

# run CUDA test
echo "  "
echo "running CUDA version, axpy.cuda:"
./axpy.cuda $N
