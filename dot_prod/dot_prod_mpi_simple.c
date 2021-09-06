/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"


/* Example routine to compute the dot-product of two vectors */
int main(int argc, char* argv[]) {

  /* declarations */
  int i, n, is, ie;
  double *a, *b, sum, mysum, alloctime, inittime, runtime;
  double stime, ftime;
  int ierr, numprocs, myid;

  /* initialize MPI */
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  /* ensure that an argument was passed in */
  if (argc < 2) {
    printf("Error: function requires one argument (vector length)\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* set n as the input argument, and ensure it's positive */
  n = atoi(argv[1]);
  if (n < 1) {
    printf("Error: vector length %i must be greater than 0\n", n);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* allocate the vectors (go ahead and put the whole thing on all procs) */
  stime = MPI_Wtime();
  a = (double *) malloc( n*sizeof(double) );
  b = (double *) malloc( n*sizeof(double) );
  ftime = MPI_Wtime();
  alloctime = ftime - stime;

  /* initialize the vector values */
  stime = MPI_Wtime();
  for (i=0; i<n; i++) {
    a[i] = (0.001 * (i + 1.0)) / n;
    b[i] = (0.001 * (n - i - 1.0)) / n;
  }
  ftime = MPI_Wtime();
  inittime = ftime - stime;

  /* start computation timer */
  stime = MPI_Wtime();

  /* determine this processor's interval */
  is = ((int) (1.0*n/numprocs))*myid;
  ie = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = n;

  /* compute dot-product */
  mysum = 0.0;
  for (i=is; i<ie; i++)  mysum += a[i]*b[i];

  /* root node collects result */
  ierr = MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* stop timer */
  ftime = MPI_Wtime();
  runtime = ftime-stime;

  /* output computed value and runtime */
  if (myid == 0) {
    printf(" vector length = %i\n",n);
    printf("   dot-product = %.16e\n",sum);
    printf("    alloc time = %.2e\n",alloctime);
    printf("     init time = %.2e\n",inittime);
    printf("      run time = %.16e\n",runtime);
  }

  /* free vectors */
  free(a);
  free(b);

  /* finalize MPI */
  ierr = MPI_Finalize();

} /* end main */

