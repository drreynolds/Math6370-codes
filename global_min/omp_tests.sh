#!/bin/bash

# build executables
make clean
make all

# run serial version
echo "  "
echo "serial version (fortran):"
./glob_min_f90.exe
echo "  "
echo "serial version (c):"
./glob_min_c.exe
echo "  "
echo "serial version (c++):"
./glob_min_cpp.exe
echo "  "

# run OpenMP versions with varying numbers of threads
for i in 1 2 3 4 5 6 7 8; do
   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP version, $i threads (fortran):"
   ./glob_min_omp_f90.exe
   echo "  "
   echo "OpenMP version, $i threads (c):"
   ./glob_min_omp_c.exe
   echo "  "
   echo "OpenMP version, $i threads (c++):"
   ./glob_min_omp_cpp.exe
   echo "  "
   echo "Auto-parallelized version, $i threads (fortran):"
   ./glob_min_auto_f90.exe
   echo "  "
   echo "Auto-parallelized version, $i threads (c):"
   ./glob_min_auto_c.exe
   echo "  "
   echo "Auto-parallelized version, $i threads (c++):"
   ./glob_min_auto_cpp.exe
done

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1


# run improved serial version
echo "  "
echo "improved serial version (fortran):"
./glob_min2_f90.exe
echo "  "
echo "improved serial version (c):"
./glob_min2_c.exe
echo "  "
echo "improved serial version (c++):"
./glob_min2_cpp.exe
echo "  "

