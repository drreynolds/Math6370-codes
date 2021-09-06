#!/bin/bash

# build executables
make clean
make serial mpi

# set number of subintervals
N=50000000
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

# run MPI versions with varying numbers of tasks
for nt in 1 2 3 4; do

   echo "  "
   echo "  "
   echo "MPI simple version, $nt tasks (fortran):"
   mpiexec -n $nt ./pi_comp_mpi_simple_f90.exe < nval
   echo "  "
   echo "MPI simple version, $nt tasks (c):"
   mpiexec -n $nt ./pi_comp_mpi_simple_c.exe < nval
   echo "  "
   echo "MPI simple version, $nt tasks (c++):"
   mpiexec -n $nt ./pi_comp_mpi_simple_cpp.exe < nval

   echo "  "
   echo "  "
   echo "MPI full version, $nt tasks (fortran):"
   mpiexec -n $nt ./pi_comp_mpi_f90.exe < nval
   echo "  "
   echo "MPI full version, $nt tasks (c):"
   mpiexec -n $nt ./pi_comp_mpi_c.exe < nval
   echo "  "
   echo "MPI full version, $nt tasks (c++):"
   mpiexec -n $nt ./pi_comp_mpi_cpp.exe < nval

done

# remove temporary stdin file
\rm nval

