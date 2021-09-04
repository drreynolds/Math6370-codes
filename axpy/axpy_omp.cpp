/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_time.h"


// Example routine to perform some simple vector linear combinations. 
int main(int argc, char* argv[]) {

  // local variables
  int n, i;
  double *x, *y, *z, a, zmax, zmaxl, stime, ftime;

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

  // allocate the vectors 
  x = new double[n];
  y = new double[n];
  z = new double[n];

  // start timer 
  stime = get_time();

  // start OpenMP parallelism 
  a = -3.0;
  zmax = -1.0e300;
  #pragma omp parallel default(shared) private(i,zmaxl) 
  {

    // output parallelism information 
#ifdef _OPENMP
    #pragma omp single
    printf(" starting OpenMP with %i processes\n", omp_get_num_threads());
#endif

    // initialize x and y 
    #pragma omp for
    for (i=0; i<n; i++) {
      x[i] = exp(2.0*(i+1)/n);
      y[i] = 1.0*(n-1)/n;
    }

    // perform linear combinations, accumulate local maximum
    zmaxl = -1.0e300;
    #pragma omp for 
    for (i=0; i<n; i++) {
      z[i] = a*x[i] + y[i];
      x[i] = y[i]/a - z[i];
      y[i] = x[i]*y[i]/n;
      zmaxl = (zmaxl > z[i]) ? zmaxl : z[i];
    }

    // combine maximum values for each thread to global max 
    #pragma omp critical(overall_max)
    { zmax = (zmax > zmaxl) ? zmax : zmaxl; }

  } // end omp parallel 

  // output maximum value in z 
  printf("  max(z) = %.16e\n",zmax);

  // stop timer 
  ftime = get_time();

  // output total time 
  printf(" runtime = %.16e\n", ftime - stime);

  // free vectors 
  delete[] x;
  delete[] y;
  delete[] z;

  return 0;
} // end main 

