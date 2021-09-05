/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */


// Inclusions 
#include <stdlib.h>
#include <stdio.h>
#include "get_time.h"



// Example routine to compute the dot-product of two vectors. 
int main(int argc, char* argv[]) {

  // declarations 
  int i, n, numprocs;
  double *a, *b, sum, alloctime, inittime, runtime;
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

  // allocate the vectors 
  stime = get_time();
  a = new double[n];
  b = new double[n];
  ftime = get_time();
  alloctime = ftime - stime;

  // initialize the vector values 
  stime = get_time();
  for (i=0; i<n; i++)    a[i] = (0.001 * (i + 1.0)) / n;
  for (i=0; i<n; i++)    b[i] = (0.001 * (n - i - 1.0)) / n;
  ftime = get_time();
  inittime = ftime - stime;

  // start timer 
  stime = get_time();

  // start OpenMP parallelism 
  sum = 0.0;
  #pragma omp parallel default(shared)
  {
    // output parallelism information 
#ifdef _OPENMP
    numprocs = omp_get_num_threads();
    #pragma omp single
    printf(" starting OpenMP with %i processes\n", numprocs);
#endif
    
    // compute dot-product 
    #pragma omp for schedule(static, n/numprocs) reduction(+:sum)
    for (i=0; i<n; i++)  sum += a[i]*b[i];

  } // end omp parallel 

  // stop timer 
  ftime = get_time();
  runtime = ftime - stime;

  // output computed value and runtime 
  printf(" vector length = %i\n",n);
  printf("   dot-product = %.16e\n",sum);
  printf("    alloc time = %.2e\n",alloctime);
  printf("     init time = %.2e\n",inittime);
  printf("      run time = %.2e\n",runtime);

  // free vectors 
  delete[] a;
  delete[] b;

  return 0;
} // end main 

