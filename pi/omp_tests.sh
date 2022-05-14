#!/bin/bash

# build executables
make clean
make all

# set number of subintervals, if not specified via the command line
if [ -n "$1" ]
then
    N=$1
else
    N=50000000
fi
echo "Running tests using N = $N"
echo "$N" > nval

# run serial versions
echo "  "
echo "serial version (fortran):"
./pi_comp_f90.exe < nval
echo "  "
echo "serial version (c):"
./pi_comp_c.exe < nval
echo "  "
echo "serial version (c++):"
./pi_comp_cpp.exe < nval

# run OpenMP versions with varying numbers of threads
for i in 1 2 3 4 5 6 7 8; do

   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP simple version, $i threads (fortran):"
   ./pi_comp_omp_simple_f90.exe < nval
   echo "  "
   echo "OpenMP simple version, $i threads (c):"
   ./pi_comp_omp_simple_c.exe < nval
   echo "  "
   echo "OpenMP simple version, $i threads (c++):"
   ./pi_comp_omp_simple_cpp.exe < nval

   echo "  "
   echo "  "
   echo "OpenMP full version, $i threads (fortran):"
   ./pi_comp_omp_f90.exe < nval
   echo "  "
   echo "OpenMP full version, $i threads (c):"
   ./pi_comp_omp_c.exe < nval
   echo "  "
   echo "OpenMP full version, $i threads (c++):"
   ./pi_comp_omp_cpp.exe < nval

   echo "  "
   echo "  "
   echo "OpenMP critical version, $i threads (fortran):"
   ./pi_comp_omp_crit_f90.exe < nval
   echo "  "
   echo "OpenMP critical version, $i threads (c):"
   ./pi_comp_omp_crit_c.exe < nval
   echo "  "
   echo "OpenMP critical version, $i threads (c++):"
   ./pi_comp_omp_crit_cpp.exe < nval

done

# remove temporary stdin file
\rm nval

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1
