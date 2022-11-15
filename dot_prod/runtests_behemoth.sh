#!/bin/bash

# set desired vector length
N=20000000

# set desired numbers of OpenMP threads
THREADS=(1 2 4 8 16 32)

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.6.00/gcc-11.2.0/lib/cmake/Kokkos

# build executables
make

# run serial test
echo "  "
echo "running serial version, dot_prod.serial:"
./dot_prod.serial $N

# run OpenMP tests
for t in "${THREADS[@]}"
do
    echo "  "
    echo "running OpenMP version, dot_prod.openmp, with $t threads:"
    OMP_NUM_THREADS=$t ./dot_prod.openmp $N

done

# run CUDA test
echo "  "
echo "running CUDA version, dot_prod.cuda:"
./dot_prod.cuda $N
