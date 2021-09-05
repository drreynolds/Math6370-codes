/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   20 January 2011 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_time.h"


/* Example routine to compute the product of an m*n matrix and an n-vector. */
int main(int argc, char* argv[]) {

  // local variables
  int m, n, i, j;
  double **A, *x, *b, norm2, runtime;
  double stime, ftime;
  FILE *FID;

  /* input the sie of the system */
  printf("We will multiply a m*n matrix by an n-vector\n");
  printf("   enter m\n");
  i = scanf("%i", &m);
  printf("   enter n\n");
  i = scanf("%i", &n);
  if ((m < 1) || (n < 1)) {
    printf(" Illegal input, m = %i and n = %i must both be >= 1\n",m,n);
    return 1;
  }

  /* allocate the matrix and vectors */
  A = (double **) malloc(m * sizeof(double *));
  for (i=0; i<m; i++)  A[i] = (double *) malloc(n * sizeof(double));
  x = (double *) malloc(n * sizeof(double));
  b = (double *) malloc(m * sizeof(double));
  
  /* initialize the matrix and vectors */
  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      A[i][j] = 1.0/(1.0 + (i-j)*(i-j));
  for (i=0; i<m; i++)  b[i] = 0.0;
  for (j=0; j<n; j++)  x[j] = 1.0;

  /* start timer */
  stime = get_time();

  /* compute matrix-vector product */
  for (j=0; j<n; j++) 
    for (i=0; i<m; i++) 
      b[i] += A[i][j]*x[j];

  /* stop timer */
  ftime = get_time();
  runtime = ftime - stime;

  /* output 2-norm of product and runtime to screen */
  norm2 = 0.0;
  for (i=0; i<m; i++)  norm2 += b[i]*b[i];
  printf("       matrix size = %i x %i\n", m, n);
  printf(" 2-norm of product = %.16e\n",sqrt(norm2));
  printf("           runtime = %.16e\n",runtime);

  /* output product to file */
  FID = fopen("b.txt","w");
  for (i=0; i<m; i++)  fprintf(FID,"%.15e\n",b[i]);
  fclose(FID);
  
  /* free matrix and vectors */
  for (i=0; i<m; i++)  free(A[i]);
  free(A);
  free(x);
  free(b);

} /* end main */

