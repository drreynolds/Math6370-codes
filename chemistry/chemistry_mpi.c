/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   20 January 2011 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "mpi.h"


/* Prototypes */
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  /* intialize MPI */
  int i, ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Init = %i\n",ierr);
    return -1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_size = %i\n",ierr);
    return -1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_rank = %i\n",ierr);
    return -1;
  }

  /* set solver input parameters */
  int maxit = 1000000;
  double lam = 1.e-2;
  double eps = 1.e-10;

  /* input the number of intervals */
  int n;
  if (myid == 0) {
    printf("Enter the total number of intervals (0 quits):\n");
    i = scanf("%i", &n);
    if (i < 1) {
      ierr = MPI_Finalize();
      return -1;
    }
  }

  /* root node outputs parallelism information to screen */
  if (myid == 0)  printf(" starting MPI with %i processes\n",numprocs);

  /* root node sends n out to other processors */
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in mpi_bcast = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }

  /* stop for illegal n */
  if (n < 1) {
    ierr = MPI_Finalize();
    return -1;
  }
  
  /* determine this processor's interval */
  int is = ((int) 1.0*n/numprocs)*myid;
  int ie = ((int) 1.0*n/numprocs)*(myid+1) - 1;
  if (myid == numprocs-1)  ie = n-1;

  /* root node allocates temperature field and solution arrays */
  double *T=NULL, *u=NULL, *v=NULL, *w=NULL;
  if (myid == 0) {
    T = malloc( n * sizeof(double) );
    u = malloc( n * sizeof(double) );
    v = malloc( n * sizeof(double) );
    w = malloc( n * sizeof(double) );

    /* set random temperature field, initial guesses at chemical densities */
    for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
    for (i=0; i<n; i++)  u[i] = 0.35;
    for (i=0; i<n; i++)  v[i] = 0.1;
    for (i=0; i<n; i++)  w[i] = 0.5;
  }
  /* other nodes allocate local temperature field and solution arrays */
  else {
    T = malloc( (ie-is+1) * sizeof(double) );
    u = malloc( (ie-is+1) * sizeof(double) );
    v = malloc( (ie-is+1) * sizeof(double) );
    w = malloc( (ie-is+1) * sizeof(double) );

    /* set initial guesses at chemical densities */
    for (i=0; i<(ie-is+1); i++)  u[i] = 0.35;
    for (i=0; i<(ie-is+1); i++)  v[i] = 0.1;
    for (i=0; i<(ie-is+1); i++)  w[i] = 0.5;
  }

  /* allocate gatherv/scatterv temporary arrays */
  int *counts = malloc( numprocs * sizeof(int) );
  int *displs = malloc( numprocs * sizeof(int) );

  /* start timer */
  double stime = MPI_Wtime();

  /* root node sends out portions of temperature field to all procs */
  int js, je;
  for (i=0; i<numprocs; i++) {
    js = floor(1.0*n/numprocs)*i;
    je = floor(1.0*n/numprocs)*(i+1) - 1;
    if (i == numprocs-1)  je = n-1;
    counts[i] = je-js+1;
    displs[i] = js;
  }
  ierr = MPI_Scatterv(T, counts, displs, MPI_DOUBLE, T, 
		      ie-is+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  /* call solver over local intervals */
  int its, it;
  double res;
  for (i=0; i<(ie-is+1); i++) {
    it = is+i;
    chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, maxit, &its, &res);
    if (res < eps) {
/*       printf("    i = %i,  its = %i\n", i, its); */
    }
    else {
      printf("    error: i = %i,  its = %i,  res = %.2e,  u = %.2e,  v = %.2e,  w = %.2e\n", i, its, res, u[i], v[i], w[i]);
      ierr = MPI_Abort(MPI_COMM_WORLD, -1);
      return -1;
    }
  }

  /* root node collects results */
  ierr = MPI_Gatherv(u, ie-is+1, MPI_DOUBLE, u, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }
  ierr = MPI_Gatherv(v, ie-is+1, MPI_DOUBLE, v, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }
  ierr = MPI_Gatherv(w, ie-is+1, MPI_DOUBLE, w, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    printf("Error in MPI_Gatherv = %i\n",ierr);
    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }

  /* stop timer */
  double ftime = MPI_Wtime();
  double runtime = ftime - stime;

  /* output solution time */
  printf(" proc %i computed %i iterations\n", myid, ie-is+1);
  printf("     runtime = %.16e\n",runtime);

  /* free temperature and solution arrays */
  free(T);
  free(u);
  free(v);
  free(w);
  free(counts);
  free(displs);

  /* finalize MPI */
  ierr = MPI_Finalize();

} /* end main */

