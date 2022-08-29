# Example codes (OpenMP)

## axpy

Simple example of the vector operation z <- a*x + y, where x, y and z are vectors and a is a scalar.

Key learning topics:
* `parallel` construct
* `for` construct
* combined `parallel for` construct

## dot_prod

Interactive demo showing parallelization of our earlier serial dot product of two vectors with OpenMP.

Key learning topics:
* `parallel` construct
* `for` construct
* `reduction` clause

## pi

Interactive demo showing OpenMP parallelization of our earlier serial code to approximate pi via numerical integration,

  pi = 4\int_{0}^{1} \frac{1}{1+x^2} dx

using the midpoint rule over n equal subintervals.

Key learning topics:
* `parallel` construct
* `for` construct
* `private` clause
* `reduction` clause

Then a "revisit" of this problem to perform the reduction manually using `private` variables and a `critical` clause for threads to contribute to the global result.

## omp_schedule

Simple OpenMP example showing the key variations of the `schedule` clause.

## global_min

Interactive demo working through an OpenMP parallelization of our global minimization example.

Key learning topics:
* `parallel` construct
* `for` construct
* `schedule` clause for load-balancing
* `private` variables

## chemistry

Interactive demo of the parallelization of our equilibrium chemical density computations using OpenMP.

Key learning topics:
* `parallel` construct
* `for` construct
* `schedule` clause for load-balancing
* orphan directives for subroutines within a parallel region

## advection

Interactive demo of the OpenMP parallelization for our 2D advection code.

Key learning topics:
* `parallel` construct
* `for` construct
* `default` clause
* `shared` clause
* `private` clause
* Dual parallelism approaches: outer loop vs inner loop

## omp_benchmarking

Example programs that help compute compiler/architecture-specific overheads associated with various OpenMP constructs.

## matvec

Example showing the OpenMP parallelization of our matrix-vector product demo code.

Key learning topics:
* Large parallel region
* The OpenMP timer
* `single` construct
* appropriate use of `for` construct on nested loops with nontrivial dependencies
