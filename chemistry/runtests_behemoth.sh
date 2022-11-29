#!/bin/bash

# set desired number of chemical bins
N=10000

# set desired numbers of OpenMP threads
THREADS=(1 2 4 8 16 32)

# configure build
cmake -DKokkos_DIR=/usr/local/kokkos-3.6.00/gcc-11.2.0/lib/cmake/Kokkos

# build executables
make

# run serial test
echo "  "
echo "running serial version, chemistry.serial:"
./chemistry.serial $N

# run OpenMP tests
for t in "${THREADS[@]}"
do
    echo "  "
    echo "running OpenMP version, chemistry.openmp, with $t threads:"
    OMP_NUM_THREADS=$t ./chemistry.openmp $N

done

# run CUDA test
echo "  "
echo "running CUDA version, chemistry.cuda:"
./chemistry.cuda $N

# run UVM test
echo "  "
echo "running UVM version, chemistry.cuda:"
./chemistry.uvm $N
