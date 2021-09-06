/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"


/* Prototypes */
inline double f(double a) { return (4.0 / (1.0 + a*a)); }


/* Example routine to compute pi using numerical integration via
      pi = 4 * int_0^1 1/(1+x^2) dx
   We use a simple midpoint rule for integration, over
   subintervals of fixed size 1/n, where n is a user-input parameter. */
int main(int argc, char* argv[]) {

  /* declarations */
  int i, n, is, ie;
  double h, x, pi, mypi, pi_true=3.14159265358979323846;
  double stime, ftime;
  int ierr, numprocs, myid;

  /* intialize MPI */
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Comm_size = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* root outputs parallelism information to screen, inputs the number of intervals */
  if (myid == 0) {
    printf(" Running with %i MPI tasks\n",numprocs);
    printf(" Enter the number of intervals (0 quits):\n");
    ierr = scanf("%i", &n);
    if (ierr != 1)  MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* root sends n out to other processors */
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  /* stop for illegal n */
  if (n < 1)  MPI_Abort(MPI_COMM_WORLD, 1);

  /* start timer */
  stime = MPI_Wtime();

  /* set subinterval width */
  h = 1.0 / n;

  /* determine this processor's interval */
  is = ((int) (1.0*n/numprocs))*myid + 1;
  ie = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = n;

  /* perform integration over n intervals */ 
  mypi = 0.0;
  for (i=is; i<=ie; i++) {
    x = h * (i - 0.5);
    mypi += h * f(x);
  }

  /* root node collects result */
  ierr = MPI_Reduce(&mypi, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Reduce = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* stop timer */
  ftime = MPI_Wtime();

  /* output computed value and error */
  if (myid == 0) {
    printf(" computed pi = %.16e\n",pi);
    printf("     true pi = %.16e\n",pi_true);
    printf("       error = %.16e\n",pi_true - pi);
    printf("     runtime = %.16e\n",ftime-stime);
  }

  /* finalize MPI */
  ierr = MPI_Finalize();

} /* end main */
