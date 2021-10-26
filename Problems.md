# Example codes (RAJA/CUDA)

## pi

Interactive demo showing RAJA/CUDA parallelization of our earlier serial code to compute pi via numerical integration,

  pi = 4\int_{0}^{1} \frac{1}{1+x^2} dx

using the midpoint rule over n equal subintervals

Key learning topics:

* `raja-template` build structure via CMake
* RAJA execution and reduction policies
* `RAJA::ReduceSum`` construct
* `RAJA::forall` loop construct
* calling device-specific functions


## dot_prod

Interactive demo showing RAJA/CUDA parallelization of our earlier serial dot-product of two vectors.

Key learning topics:

* `cudaMalloc` and `cudaFree` routines for memory allocation/deallocation on the device


## axpy

Simple example that computes the linear combination a*x + y, where x and y are vectors, and a is a scalar.

Key learning topics:

* `RAJA::ReduceMax` reduction construct


## global_min

Interactive demo working through a RAJA/CUDA parallelization of our previous lab on performing steepest-descent minimization repeatedly using different initial conditions.

Key learning topics:

*


## chemistry

Interactive demo of the RAJA/CUDA paralleization of our equilibrium chemical density computations.

Key learning topics:
*


## advection

Interactive demo of the RAJA/CUDA parallelization of our 2D advection code.

Key learning topics:
*
