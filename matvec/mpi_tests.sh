#!/bin/bash

# build executables
make realclean
make all

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

# run MPI versions with varying numbers of tasks
for nt in 1 2 3 4; do

   echo "  "
   echo "  "
   echo "MPI version 1, $nt tasks (fortran):"
   mpiexec -n $nt ./matvec_mpi1_f90.exe < sizes
   echo "  "
   echo "MPI version 1, $nt tasks (c):"
   mpiexec -n $nt ./matvec_mpi1_c.exe < sizes
   echo "  "
   echo "MPI version 1, $nt tasks (c++):"
   mpiexec -n $nt ./matvec_mpi1_cpp.exe < sizes

   echo "  "
   echo "  "
   echo "MPI version 2, $nt tasks (fortran):"
   mpiexec -n $nt ./matvec_mpi2_f90.exe < sizes
   echo "  "
   echo "MPI version 2, $nt tasks (c):"
   mpiexec -n $nt ./matvec_mpi2_c.exe < sizes
   echo "  "
   echo "MPI version 2, $nt tasks (c++):"
   mpiexec -n $nt ./matvec_mpi2_cpp.exe < sizes

   echo "  "
   echo "  "
   echo "MPI version 3, $nt tasks (fortran):"
   mpiexec -n $nt ./matvec_mpi3_f90.exe < sizes
   echo "  "
   echo "MPI version 3, $nt tasks (c):"
   mpiexec -n $nt ./matvec_mpi3_c.exe < sizes
   echo "  "
   echo "MPI version 3, $nt tasks (c++):"
   mpiexec -n $nt ./matvec_mpi3_cpp.exe < sizes

done

# remove temporary stdin file
\rm sizes
