/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


/* Example routine to compute the product of an m*n matrix and an n-vector. */
int main(int argc, char* argv[]) {

  /* local variables */
  int m, n, i, j, js, je, ierr, numprocs, myid;
  double **A, *x, *b, *myb, runtime, norm2;
  double stime, ftime;
  FILE *FID;

  /* intialize MPI */
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_size = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* input the sie of the system */
  if (myid == 0) {
    printf("We will multiply a m*n matrix by an n-vector\n");
    printf("   enter m\n");
    ierr = scanf("%i", &m);
    printf("   enter n\n");
    ierr = scanf("%i", &n);
  }

  /* root node sends m and n out to other processors */
  ierr = MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if ((m < 1) || (n < 1)) {
    if (myid == 0)
      fprintf(stderr," Illegal input, m = %i and n = %i must both be >= 1\n",m,n);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* root node outputs parallelism information to screen */
  if (myid == 0)
    printf(" starting MPI with %i processes\n",numprocs);

  /* determine this processor's interval */
  js = ((int) (1.0*n/numprocs))*myid;
  je = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  je = n;

  /* allocate the matrix and vectors */
  /* (store entire system on every proc -- wasteful) */
  A = (double **) malloc( m * sizeof(double *));
  for (i=0; i<m; i++)  A[i] = (double *) malloc((je-js) * sizeof(double));
  x = (double *) malloc((je-js) * sizeof(double));
  b = (double *) malloc(m * sizeof(double));
  myb = (double *) malloc(m * sizeof(double));

  /* initialize the matrix and vectors */
  for (i=0; i<m; i++)
    for (j=js; j<je; j++)
      A[i][j-js] = 1.0/(1.0 + (i-j)*(i-j));
  for (j=0; j<m; j++)    b[j] = 0.0;
  for (j=js; j<je; j++)  x[j-js] = 1.0;

  /* start timer */
  stime = MPI_Wtime();

  /* compute matrix-vector product */
  for (i=0; i<m; i++)  myb[i] = 0.0;
  for (i=0; i<m; i++)
    for (j=js; j<je; j++)
      myb[i] += A[i][j-js]*x[j-js];

  /* root node collects result */
  ierr = MPI_Reduce(myb, b, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    fprintf(stderr," error in MPI_Reduce = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* stop timer */
  ftime = MPI_Wtime();

  /* output 2-norm of product and runtime to screen */
  if (myid == 0) {
    norm2 = 0.0;
    for (i=0; i<m; i++)  norm2 += b[i]*b[i];
    printf("       matrix size = %i x %i\n", m, n);
    printf(" 2-norm of product = %.16e\n",sqrt(norm2));
    printf("           runtime = %.16e\n",ftime-stime);

    /* output product to file */
    FID = fopen("b_mpi3.txt","w");
    for (i=0; i<m; i++)  fprintf(FID,"%.15e\n",b[i]);
    fclose(FID);
  }

  /* free matrix and vectors */
  for (i=0; i<m; i++)  free(A[i]);
  free(A);
  free(x);
  free(b);
  free(myb);

  /* finalize MPI */
  ierr = MPI_Finalize();

} /* end main */
