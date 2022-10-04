# Example codes (Kokkos)

## pi

Interactive demo showing Kokkos parallelization of our earlier serial code to compute pi via numerical integration,
$$\pi = 4\int_{0}^{1} \frac{1}{1+x^2} dx$$
using the midpoint rule over n equal subintervals

Key learning topics:

* `CMake` build system, including multiple targets with preprocessor directives
* `Kokkos::initialize`
* Kokkos execution space and reduction policies
* `Kokkos::Timer`
* `Kokkos::parallel_reduce` construct


## dot_prod

Interactive demo showing Kokkos parallelization of our earlier serial dot-product of two vectors.

Key learning topics:

* Kokkos memory space and 1D views for portable host/device array allocation/deallocation
* `Kokkos::parallel_for` construct
* `Kokkos::fence()` to ensure completion of asynchronous execution


## axpy

Simple example that computes the linear combination a*x + y, where x and y are vectors, and a is a scalar.

Key learning topics:

* Named Kokkos parallel loops (beneficial for profiling/debugging)
* "Maximum" reduction


## global_min

Interactive demo working through a Kokkos parallelization of our previous lab on performing steepest-descent minimization repeatedly using different initial conditions.

Key learning topics:

* Definition of external functions that are callable from either CPU or GPU code.
* Custom reducer for `Kokkos::parallel_reduce` -- allows us to identify (x,y) leading to global minimum on GPU, with result automatically transferred to CPU.


## chemistry

Interactive demo of the Kokkos parallelization of our equilibrium chemical density computations.

Key learning topics:

* `cudaMemcpy` for host->device transfer of temperature array, and then device->host transfer of chemical equilibrium results.


## advection

Interactive demo of the Kokkos parallelization of our 2D advection code.

Key learning topics:

* `cudaMemcpyAsync` for transfers of device solution data to host for output to disk.
* RAJA "view"
