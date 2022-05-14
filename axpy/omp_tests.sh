#!/bin/bash

# build executables
make realclean
make all

# set vector length
if [ -n "$1" ]
then
    N=$1
else
    N=40000000
fi
echo "Running tests using N = $N"


# run serial versions
echo "  "
echo "  "
echo "serial version (fortran):"
./axpy_f90.exe $N
echo "  "
echo "serial version (c):"
./axpy_c.exe $N
echo "  "
echo "serial version (c++):"
./axpy_cpp.exe $N

# run auto-parallelized and OpenMP versions with varying numbers of threads
for i in 1 2 3 4 5 6 7 8; do
   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP version, $i threads (fortran):"
   ./axpy_omp_f90.exe $N
   echo "  "
   echo "OpenMP version, $i threads (c):"
   ./axpy_omp_c.exe $N
   echo "  "
   echo "OpenMP version, $i threads (c++):"
   ./axpy_omp_cpp.exe $N

   if [ $i -eq 4 ]
   then
       echo "  "
       echo "  "
       echo "Auto-parallelized version, $i threads (fortran):"
       ./axpy_auto_f90.exe $N
       echo "  "
       echo "Auto-parallelized version, $i threads (c):"
       ./axpy_auto_c.exe $N
       echo "  "
       echo "Auto-parallelized version, $i threads (c++):"
       ./axpy_auto_cpp.exe $N
   fi
done

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1
