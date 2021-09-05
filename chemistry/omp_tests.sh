#!/bin/bash

# build executables
make clean
make all

# set number of subintervals
N=2000
echo "$N" > nval

# run OpenMP versions with varying numbers of threads
for i in 1 2 3 4; do

   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP version, $i threads (fortran):"
   ./chemistry_omp_f90.exe < nval
   echo "  "
   echo "OpenMP version, $i threads (c):"
   ./chemistry_omp_c.exe < nval
   echo "  "
   echo "OpenMP version, $i threads (c++):"
   ./chemistry_omp_cpp.exe < nval

   echo "  "
   echo "  "
   echo "OpenMP orphan version, $i threads (fortran):"
   ./chemistry_omp2_f90.exe < nval
   echo "  "
   echo "OpenMP orphan version, $i threads (c):"
   ./chemistry_omp2_c.exe < nval
   echo "  "
   echo "OpenMP orphan version, $i threads (c++):"
   ./chemistry_omp2_cpp.exe < nval

done

# remove temporary stdin file
\rm nval

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1

