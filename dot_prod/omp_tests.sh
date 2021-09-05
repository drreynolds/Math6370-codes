#!/bin/bash

# build executables
make clean
make all

# set vector length
N=50000000

# run serial versions
echo "  "
echo "  "
echo "serial version (fortran):"
./dot_prod_f90.exe $N
echo "  "
echo "serial version (c):"
./dot_prod_c.exe $N
echo "  "
echo "serial version (c++):"
./dot_prod_cpp.exe $N

# run OpenMP versions with varying numbers of threads
for i in 1 2 3 4 5 6 7 8; do
   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP version, $i threads (fortran):"
   ./dot_prod_omp_f90.exe $N
   echo "  "
   echo "OpenMP version, $i threads (c):"
   ./dot_prod_omp_c.exe $N
   echo "  "
   echo "OpenMP version, $i threads (c++):"
   ./dot_prod_omp_cpp.exe $N

   echo "  "
   echo "  "
   echo "OpenMP version 2, $i threads (fortran):"
   ./dot_prod_omp2_f90.exe $N
   echo "  "
   echo "OpenMP version 2, $i threads (c):"
   ./dot_prod_omp2_c.exe $N
   echo "  "
   echo "OpenMP version 2, $i threads (c++):"
   ./dot_prod_omp2_cpp.exe $N
done

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1
