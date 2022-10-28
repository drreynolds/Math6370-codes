# HYPRE MPI example code

## Examples

This directory includes 4 separate HYPRE examples, all of which use HYPRE's structure-grid interface.  Each example sets up and solves a scalar-valued Poisson-like problem, with homogeneous Dirichlet boundary conditions,
$$(I + L)v = w,$$
where $L$ is a standard two-dimensional, 5-point Laplace stencil operator (with potentially anisotropic diffusion coefficients), and $w$ is a smooth right-hand side vector.

* `hypre_test.cpp`: the Laplace operator is discretized using a simple second-order centered difference approximation on a cell-centered finite-volume grid.  The problem is solved using HYPRE's PFMG geometric multigrid solver.

* `hypre_test2.cpp`: the Laplace operator is discretized using a simple second-order centered difference approximation on a cell-centered finite-volume grid.  The problem is solved using HYPRE's BiCGStab linear solver, that is preconditioned using HYPRE's PFMG geometric multigrid solver.  We use the input file parameter `PCGmaxit` for the number of outer BiCGStab iterations.

* `hypre_test_fd.cpp`: the Laplace operator is discretized using a simple second-order centered difference approximation on a nodal finite-difference grid.  The problem is solved using HYPRE's PFMG geometric multigrid solver.

  Note that although we might naturally choose to require a grid size of Nx=$$2^k+1$$ for geometric multigrid solvers on finite-difference grids, HYPRE's PFMG geometric-multigrid solver requires that Nx=$$2^k$$.

* `hypre_test2_fd.cpp`: the Laplace operator is discretized using a simple second-order centered difference approximation on a cell-centered finite-volume grid.  The problem is solved using HYPRE's PCG linear solver, that is preconditioned using HYPRE's PFMG geometric multigrid solver.

  As with `hypre_test_fd.cpp`, we must use a grid with sizes Nx=$$2^k$$.


## Building

This example requires compilation against a correctly-installed version of the [HYPRE linear solver library](https://github.com/hypre-space/hypre).  We are uncertain of the precise *minimum* version of HYPRE that is required for our usage; however, any installation of version 2.16.0 or higher should work.

Since this may be installed in any location on a given system, the current `Makefile` relies on environment variables to build this example:

* `HYPRE_ROOT`: the top-level directory for the HYPRE installation to use.  This must have the subdirectories `include` and `lib`, containing the HYPRE header and library files, respectively.

* `MPICXX`: the MPI C++ compiler wrapper that corresponds with the same MPI library as was used to compile the HYPRE library.

We also recommend that if the `mpiexec` corresponding with the `MPICXX` above is *not* the first one in your `PATH`, that you also set an environment variable

* `MPIEXEC`: the MPI execution launch script that corresponds with `MPICXX` above.

[Daniel R. Reynolds](https://github.com/drreynolds)
[Mathematics @ SMU](https://www.smu.edu/math)
