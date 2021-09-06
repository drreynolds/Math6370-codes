# Example codes (MPI)

## simple

Simple example showing the basic structure of an MPI program.

## dot_prod

Interactive demo showing parallelization of our earlier serial dot product of two vectors with MPI.

Versions and key learning topics:
* `dot_prod_mpi_simple`:
  * full vectors are stored on all MPI ranks
  * uses `MPI_Reduce` to combine contributions together
  * switches to MPI timers
* `dot_prod_mpi`: adds simple error handling to `dot_prod_mpi_simple`
* `dot_prod_mpi_fancy`: converts `dot_prod_mpi` to only allocate rank-specific portions of input vectors.
* `dot_prod_mpi_noreduce`: converts `dot_prod_mpi_fancy` to perform the reduction manually (instead of calling `MPI_Reduce`).

## pi

Interactive demo showing MPI parallelization of our earlier serial code to approximate pi via numerical integration,

  pi = 4\int_{0}^{1} \frac{1}{1+x^2} dx

using the midpoint rule over n equal subintervals.

Versions and key learning topics:
* `pi_comp_mpi_simple`:
  * root process input n from stdin and calls `MPI_Bcast` to send value to other ranks.
  * uses `MPI_Reduce` to combine contributions together
* `pi_comp_mpi`: adds simple error handling to `pi_comp_mpi_simple`

## matvec

In-class demo, MPI parallelizing serial matrix-vector product demo code.

Versions and key learning topics:
* `matvec_mpi1`:
  * Root process inputs m and n from stdin and calls `MPI_Bcast` to send values to other ranks.
  * Stores full matrix on every rank; only compute your own rows.
  * Root process collects overall result using `MPI_Gather`.
  * Simple error handling in case calls fail.
* `matvec_mpi2`: converts `matvec_mpi1` to only store a "block row decomposition" of the matrix, where each rank only stores its local matrix rows.
* `matvec_mpi3`: converts `matvec_mpi1` to only store a "block column decomposition" of the matrix, where each rank only stores its local matrix columns.

## laplace

Linear residual calculation demo for a 2D Laplace operator on a regular finite difference grid.

`linresid2D_mpi_sync`:
* Point-to-point communication is based on `MPI_Send` and `MPI_Recv`:
* Requires an 8-phase system for 2D communication to neighbors to avoid deadlock due to an incorrect ordering of messages:
  * X-directional communication:
    * Evens send East, while odds receive West
    * Odds send West, evens receive East
    * Odds send East, evens send West
    * Evens send West, odds receive East
  * Y-directional communication:
    * Evens send North, odds receive South
    * Odds send South, evens receive North
    * Odds send North, evens receive South
    * Evens send South, odds receive North
  * Once all data is available, perform the residual computation.

`linresid2D_mpi_async`:
* Point-to-point communication is based on `MPI_Isend` and `MPI_Irecv`:
* Instead of our complicated 8-phase system for 2D communication to neighbors, we may
communicate boundary information to neighbors in a simpler manner. Each process will:
  * Allocate receive buffers for neighbor information,
  * Open receive channels to fill these buffers, using `MPI_Irecv`,
  * Fill send buffers to send to neighbors,
  * Send information to waiting neighbor processes, using `MPI_Isend`,
  * Wait for the information to arrive, using `MPI_Wait`,
  * Perform the residual computation.
* We no longer have to worry about deadlock due to an incorrect ordering of messages.
* We could even do some computations while we wait for the boundary data to arrive.
* This overlapping of communication and computation is desirable, since codes can
operate at speeds much faster than the communication network.
* We must be careful with the tags, since each process may receive multiple messages
simultaneously, and we need to ensure that they end up in the right places.

## HYPRE

Example showing how to use HYPRE to perform a scalable geometric multigrid linear solve for a linear system arising from spatial discretization of a 2D Poisson equation.

## advection

Interactive demo of the MPI parallelization for our 2D advection code.

Versions and key learning topics:
* `*_mpi`:
  * Periodic Cartesian communicator
  * Asynchronous communication, using `MPI_Waitall`
* `*_1sided`: MPI3+ version using one-sided communication: `MPI_Put` and "windows".

## chemistry

Interactive MPI parallelization of our equilibrium chemical density computations.

* `chemistry_mpi` --
  * Static/even work distribution among MPI ranks
  * Root process assigns temperatures via `MPI_Scatterv` and obtains results using `MPI_Gatherv`
* `chemistry_mpi2` --
  * Manager/worker distribution among MPI ranks (one bin at a time)
* `chemistry_mpi3` --
  * Manager/worker distribution among MPI ranks (10 bins at a time)
* `chemistry_hybrid.cpp` -- hybrid MPI+OpenMP parallelization of `chemistry_mpi` code
  * Initialize MPI in threadsafe manner via `MPI_Init_thread`
  * Exploration of hybrid parallelism performance: pure MPI, pure OpenMP, hybrid MPI+OpenMP.
