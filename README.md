# Math6370-codes

Codes for in-class collaboration for the course: Parallel Scientific Computing (MATH 4370/6370) at [Southern Methodist University](https://www.smu.edu).

These codes require working C, C++11, and Fortran-90 compilers that support OpenMP, a working MPI installation (preferably with the `mpicc`, `mpicxx` and `mpifort` wrapper scripts and `mpiexec` launcher), and the "make" program, and are designed to be compiled/run at the Linux/OSX command line.  For the `raja` branch, we additionally require that [RAJA](https://github.com/LLNL/RAJA) v0.12.1 or higher is installed, with the [CUDA](https://developer.nvidia.com/cuda-downloads) backend (this introduces additional dependencies).  Similarly, for the `kokkos` branch we require that [Kokkos](https://github.com/kokkos/kokkos) v3.6.01 or higher is installed.  These codes have been tested with the [CUDA](https://developer.nvidia.com/cuda-downloads) backend and `cuda_lambda` option enabled; however, we hope that the codes will also work with AMD or Intel GPUs and the HIP or SYCL backends, respectively.  

On most Linux systems, MPI, RAJA and Kokkos (and their dependencies) may be easily installed using [Spack](https://github.com/spack/spack).  

Note: we use branches to focus on specific code types:
* main: serial
* openmp: shared-memory parallelization via OpenMP
* mpi: distributed-memory parallelization via MPI
* raja: accelerator-based parallelism with RAJA (C++ only)
* kokkos: accelerator-based parallelism with Kokkos (C++ only)

In the `kokkos` branch, we provide `runtests_m2.sh` scripts in each
example folder.  These assume that you are running on the ManeframeII
cluster at SMU, that you have already acquired an interactive session
on a GPU node,
```
srun -p development -c 4 --mem=16G --gres=gpu:volta:1 --pty $SHELL
```
and that you have already configured your environment for using Kokkos,
```
module load spack gcc-9.2
. /hpc/spack/share/spack/setup-env.sh
spack load kokkos/qu45u5v
```

[Daniel R. Reynolds](https://github.com/drreynolds)  
[Mathematics @ SMU](https://www.smu.edu/math)
