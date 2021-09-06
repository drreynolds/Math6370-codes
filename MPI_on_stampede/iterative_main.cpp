/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   28 February 2013 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


// Prototypes
double linresid(double *, double *, double *, double *, 
		double *, double *, int, int, MPI_Comm);
double jacobi_solve(double *, double *, double *, double *, double *, 
		    double *, int, int, double, int, int *, MPI_Comm);


/* We set up and solve the linear system 
             (I + gamma*L)u = r,
   where L is a standard 1D Laplace operator, r is a given 
   right-hand side, and u is the solution, using a parallelized 
   Jacobi iterative solver.
     
   Requires three input arguments: gamma, delta, and n 
   (the size the global domain). */
int main(int argc, char* argv[]) {

  // declarations
  int ierr, nprocs, my_id, global_N, local_N, k, maxiter, iters;
  double gamma, delta, err2norm, stime, ftime;
  double *u, *r, *a, *b, *c, *res;
  FILE *FID;
  MPI_Comm comm;

  // intialize MPI
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Init = %i\n",ierr);
    return 1;
  }
  comm = MPI_COMM_WORLD;
  ierr = MPI_Comm_size(comm, &nprocs);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_size = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Comm_rank(comm, &my_id);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // root node gets problem information from input file
  if (my_id == 0) {
    FID = fopen("input_iterative.txt","r");
    fscanf(FID," &inputs\n");
    fscanf(FID,"  gamma = %lf,\n", &gamma);
    fscanf(FID,"  delta = %lf,\n", &delta);
    fscanf(FID,"  global_N = %i,\n", &global_N);
    fclose(FID);
  }

  // root node broadcasts information to other procs
  ierr = MPI_Bcast(&gamma, 1, MPI_DOUBLE, 0, comm);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Bcast(&delta, 1, MPI_DOUBLE, 0, comm);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Bcast(&global_N, 1, MPI_INT, 0, comm);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // set local mesh sizes (last process takes up remainder)
  local_N = global_N/nprocs;
  if (my_id == nprocs-1) 
    if (local_N*nprocs != global_N) 
      local_N = global_N - local_N*(nprocs-1);

  // root node outputs some information to screen
  if (my_id == 0) {
    printf("iterative test with %i processors\n",nprocs);
    printf("    gamma = %g\n",gamma);
    printf("    linear solver tolerance delta = %g\n",delta);
    printf("    global problem size N = %i\n",global_N);
    printf("    local problem sizes n = %i\n",local_N);
  }

  // Allocate vector memory
  u   = new double[local_N];
  r   = new double[local_N];
  a   = new double[local_N];
  b   = new double[local_N];
  c   = new double[local_N];
  res = new double[local_N];

  // Set up matrix arrays, right-hand side, and initial solution guess
  for (k=0; k<local_N; k++) {
    u[k] =  0.0;
    r[k] =  1.0;
    a[k] = -gamma;
    b[k] =  1.0+gamma*2.0;
    c[k] = -gamma;
  }
  
  // Adjust a, c arrays if we are at either end of domain
  if (my_id == 0)         a[0] = 0.0;
  if (my_id == nprocs-1)  c[local_N-1] = 0.0;

  // check linear residual
  err2norm = linresid(a, b, c, u, r, res, local_N, global_N, comm);
  if (my_id == 0)  
    printf(" initial residual: ||T*u-r||_2 = %g\n", err2norm);

  // Wait until all procs have caught up
  ierr = MPI_Barrier(comm);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Barrier = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // Solve system, get timing information
  maxiter = 10000;
  stime = MPI_Wtime();
  err2norm = jacobi_solve(a, b, c, u, r, res, local_N, global_N, 
			  delta, maxiter, &iters, comm);
  ftime = MPI_Wtime();
  if (my_id == 0) {
    printf(" converged in %i iterations at delta = %g\n", iters, delta);
    printf(" reported residual: ||T*u-r||_2 = %g\n", err2norm);
    printf(" solution time = %g seconds\n", ftime-stime);
  }
  
  // check linear residual
  err2norm = linresid(a, b, c, u, r, res, local_N, global_N, comm);
  if (my_id == 0)  
    printf(" actual final residual: ||T*u-r||_2 = %g\n", err2norm);

  // Free matrix/solver memory
  delete[] u;
  delete[] r;
  delete[] a;
  delete[] b;
  delete[] c;
  delete[] res;

  // finalize MPI
  MPI_Finalize();

  return 0;

} // end main

