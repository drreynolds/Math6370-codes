#!/bin/bash

# build executables
make realclean
make serial openmp

# set number of subintervals
M=20000
N=10000
echo "$M" > sizes
echo "$N" >> sizes

# run serial versions
echo "  "
echo "serial, row-oriented version (fortran):"
./matvec_row_f90.exe < sizes
echo "  "
echo "serial, column-oriented version (fortran):"
./matvec_col_f90.exe < sizes
echo "  "
echo "serial, row-oriented version (c):"
./matvec_row_c.exe < sizes
echo "  "
echo "serial, column-oriented version (c):"
./matvec_col_c.exe < sizes
echo "  "
echo "serial, row-oriented version (c++):"
./matvec_row_cpp.exe < sizes
echo "  "
echo "serial, column-oriented version (c++):"
./matvec_col_cpp.exe < sizes

# run OpenMP versions with varying numbers of threads
for i in 1 2 3 4; do

   export OMP_NUM_THREADS=$i
   echo "  "
   echo "  "
   echo "OpenMP version (column-oriented), $i threads (fortran):"
   ./matvec_omp_f90.exe < sizes
   echo "  "
   echo "OpenMP version (row-oriented), $i threads (c):"
   ./matvec_omp_c.exe < sizes
   echo "  "
   echo "OpenMP version (row-oriented), $i threads (c++):"
   ./matvec_omp_cpp.exe < sizes

done

# remove temporary stdin file
\rm sizes

# reset OMP_NUM_THREADS
export OMP_NUM_THREADS=1
