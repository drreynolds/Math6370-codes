/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_time.h"


/* Example routine to perform some simple vector linear combinations. */
int main(int argc, char* argv[]) {

  /* set the vector length */
  int n, i;
  double *x, *y, *z, a, zmax, runtime;
  double stime, ftime;
  
   // ensure that an argument was passed in
  if (argc < 2) {
    printf("Error: function requires one argument (vector length)\n");
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  n = atoi(argv[1]);
  if (n < 1) {
    printf("Error: vector length %i must be greater than 0\n", n);
    return 1;
  }

  /* allocate the vectors */
  x = malloc(n*sizeof(double));
  y = malloc(n*sizeof(double));
  z = malloc(n*sizeof(double));

  /* start timer */
  stime = get_time();

  /* initialize a, x and y */
  a = -3.0;
  for (i=0; i<n; i++) {
    x[i] = exp(2.0*(i+1)/n);
    y[i] = 1.0*(n-1)/n;
  }

  /* perform linear combinations */
  for (i=0; i<n; i++)  z[i] = a*x[i] + y[i];
  for (i=0; i<n; i++)  x[i] = y[i]/a - z[i];
  for (i=0; i<n; i++)  y[i] = x[i]*y[i]/n;

  /* output maximum value in z */
  zmax = z[0];
  for (i=0; i<n; i++) zmax = (zmax > z[i]) ? zmax : z[i];
  printf("  max(z) = %.16e\n",zmax);

  /* stop timer */
  ftime = get_time();
  runtime = ftime - stime;

  /* output total time */
  printf(" runtime = %.16e\n",runtime);

  /* free vectors */
  free(x);
  free(y);
  free(z);

  return 0;
} /* end main */

