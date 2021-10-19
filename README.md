# Math6370-codes

Codes for in-class collaboration for the course: Parallel Scientific Computing (MATH 4370/6370) at Southern Methodist University.

These codes require working C, C++11, and Fortran-90 compilers that support OpenMP, a working MPI installation (preferably with the `mpicc`, `mpicxx` and `mpifort` wrapper scripts and `mpiexec` launcher), and the "make" program, and are designed to be compiled/run at the Linux/OSX command line.

Note: we use branches to focus on specific code types:
* main: serial
* openmp: shared-memory parallelization via OpenMP
* mpi: distributed-memory parallelization via MPI
* raja: accelerator-based parallelism with RAJA

Daniel R. Reynolds  
Mathematics @ SMU
