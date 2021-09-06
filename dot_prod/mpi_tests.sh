#!/bin/bash

# build executables
make clean
make serial mpi

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

# run MPI versions with varying numbers of threads
for nt in 1 2 3 4; do
   echo "  "
   echo "  "
   echo "MPI version, $nt procs (fortran):"
   mpiexec -n $nt ./dot_prod_mpi_f90.exe $N
   echo "  "
   echo "MPI version, $nt procs (c):"
   mpiexec -n $nt ./dot_prod_mpi_c.exe $N
   echo "  "
   echo "MPI version, $nt procs (c++):"
   mpiexec -n $nt ./dot_prod_mpi_cpp.exe $N

   echo "  "
   echo "  "
   echo "MPI version (fancy), $nt procs (fortran):"
   mpiexec -n $nt ./dot_prod_mpi_fancy_f90.exe $N
   echo "  "
   echo "MPI version (fancy), $nt procs (c):"
   mpiexec -n $nt ./dot_prod_mpi_fancy_c.exe $N
   echo "  "
   echo "MPI version (fancy), $nt procs (c++):"
   mpiexec -n $nt ./dot_prod_mpi_fancy_cpp.exe $N
done

