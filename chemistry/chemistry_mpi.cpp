/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   4 January 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


// Prototypes
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // declarations
  int maxit, n, i, its, it, is, ie, js, je;
  double lam, eps, *T, *u, *v, *w, res, runtime;
  double stime, ftime;
  int ierr, numprocs, myid;
  int *counts, *displs;

  // intialize MPI
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_size = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_rank = %i\n",ierr);
    return 1;
  }

  // set solver input parameters
  maxit = 1000000;
  lam   = 1.e-2;
  eps   = 1.e-10;

  // root node outputs parallelism information to screen
  if (myid == 0)  printf("Starting MPI with %i processes\n",numprocs);

  // root proc inputs the number of intervals
  if (myid == 0) {
    printf("Enter the total number of intervals (0 quits):\n");
    i = scanf("%i", &n);
    if (i < 1) {
      ierr = MPI_Finalize();
      return 1;
    }
  }

  // root node sends n to other processors
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Bcast = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // stop for illegal n
  if (n < 1) {
    ierr = MPI_Finalize();
    return 1;
  }
  
  // determine this processor's interval
  is = ((int) 1.0*n/numprocs)*myid;
  ie = ((int) 1.0*n/numprocs)*(myid+1) - 1;
  if (myid == numprocs-1)  ie = n-1;

  // root node allocates temperature field and solution arrays
  if (myid == 0) {
    T = new double[n];
    u = new double[n];
    v = new double[n];
    w = new double[n];

    // set random temperature field, initial guesses at chemical densities
    for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
    for (i=0; i<n; i++)  u[i] = 0.35;
    for (i=0; i<n; i++)  v[i] = 0.1;
    for (i=0; i<n; i++)  w[i] = 0.5;
  }
  // other nodes allocate local temperature field and solution arrays
  else {
    T = new double[ie-is+1];
    u = new double[ie-is+1];
    v = new double[ie-is+1];
    w = new double[ie-is+1];

    // set initial guesses at chemical densities
    for (i=0; i<(ie-is+1); i++)  u[i] = 0.35;
    for (i=0; i<(ie-is+1); i++)  v[i] = 0.1;
    for (i=0; i<(ie-is+1); i++)  w[i] = 0.5;
  }

  // allocate gatherv/scatterv temporary arrays
  counts = new int[numprocs];
  displs = new int[numprocs];

  // start timer
  stime = MPI_Wtime();

  // root node sends out portions of temperature field to all procs
  for (i=0; i<numprocs; i++) {
    js = floor(1.0*n/numprocs)*i;
    je = floor(1.0*n/numprocs)*(i+1) - 1;
    if (i == numprocs-1)  je = n-1;
    counts[i] = je-js+1;
    displs[i] = js;
  }
  ierr = MPI_Scatterv(T, counts, displs, MPI_DOUBLE, T, 
		      ie-is+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  // everyone calls solver over local intervals
  for (i=0; i<(ie-is+1); i++) {
    it = is+i;
    chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, maxit, &its, &res);
    if (res < eps) {
      // printf("    i = %i,  its = %i\n", i, its);
    }
    else {
      printf("    error: i = %i,  its = %i,  res = %.2e,  u = %.2e,  v = %.2e,  w = %.2e\n", 
	     i, its, res, u[i], v[i], w[i]);
      ierr = MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
  }

  // root node collects results
  ierr = MPI_Gatherv(u, ie-is+1, MPI_DOUBLE, u, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  ierr = MPI_Gatherv(v, ie-is+1, MPI_DOUBLE, v, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  ierr = MPI_Gatherv(w, ie-is+1, MPI_DOUBLE, w, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // stop timer
  ftime = MPI_Wtime();
  runtime = ftime - stime;

  // output solution time on each proc
  printf(" proc %i computed %i iterations, runtime = %.16e\n", 
	 myid, ie-is+1, runtime);

  // free temperature and solution arrays
  delete[] T;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] counts;
  delete[] displs;

  // finalize MPI
  ierr = MPI_Finalize();

} // end main
