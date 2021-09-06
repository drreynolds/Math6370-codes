/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   28 February 2013 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


// Function Prototypes
double linresid(double *, double *, double *, double *, 
		double *, double *, int, int, MPI_Comm);


/* Description: 
     Solves the linear system 
             a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k),
     using a parallelized Jacobi iterative solver.  
     Returns the final WRMS norm of the linear residual.  */
double jacobi_solve(double *a, double *b, double *c, double *u, 
		    double *r, double *res, int local_N, int global_N, 
		    double delta, int maxiter, int *iters, MPI_Comm comm) {

  // declarations
  int ierr, i, its;
  double resid2;

  // compute initial linear residual
  resid2 = linresid(a, b, c, u, r, res, local_N, global_N, comm);

  // iterate until resid2 < delta (or maxiter iterations)
  for (its=0; its<=maxiter; its++) {

    // check for convergence
    if (resid2 < delta)  break;

    // update u = u - diag(b)^{-1}*res
    for (i=0; i<local_N; i++)   u[i] = u[i] - res[i]/b[i];

    // compute linear residual
    resid2 = linresid(a, b, c, u, r, res, local_N, global_N, comm);
  } // for its

  if (its == maxiter)
    fprintf(stderr," jacobi_solve warning: reached maximum iteration limit!\n");

  *iters = its;
  return resid2;

} // end jacobi_solve
