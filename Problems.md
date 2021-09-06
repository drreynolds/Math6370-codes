# Example codes (serial)

## dot_prod

Simple example that computes the dot product of two vectors.

## pi

Example that approximates pi via numerical integration,

  pi = 4\int_{0}^{1} \frac{1}{1+x^2} dx

using the midpoint rule over n equal subintervals

## chemistry

Example that computes equilibrium chemical densities at multiple spatial locations, given a random background temperature field, using a simple damped fixed-point iteration as the nonlinear solver.

Key learning topics:
* Calling external functions
* Call by reference vs call by value
* Standard math library
* Control structures

## advection

This example sets up and evolves the 2D first order wave equations in time, using a regular staggered spatial discretization and an explicit "leapfrog" time discretization.

Key learning topics:
* File input/output
* Header files, macros
* "Flattened" 2D -> 1D array data structures
* Nontrivial loop structures
* "() ? () : ()" conditional operators
* Python postprocessing/visualization of simulation output data

## memory_leak

This example demonstrates how we can use the Linux program "valgrind" to help identify the location of memory leaks in a program.

Key learning topics:
* Memory leaks
* Valgrind

## matvec

This example performs a simple matrix-vector product, but where we access the matrix data in different orders, illustrating that knowledge of memory layout is crucial for efficient programs.

Key learning topic: matching loop structure to data layout
