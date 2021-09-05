/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "get_time.h"



/* Example routine to compute the dot-product of two vectors. */
int main(int argc, char* argv[]) {

  /* declarations */
  int i, n, numprocs;
  double *a, *b, sum, alloctime, inittime, runtime;
  double stime, ftime;

  /* ensure that an argument was passed in */
  if (argc < 2) {
    printf("Error: function requires one argument (vector length)\n");
    return 1;
  }
    
  /* set n as the input argument, and ensure it's positive */
  n = atoi(argv[1]);
  if (n < 1) {
    printf("Error: vector length %i must be greater than 0\n", n);
    return 1;
  }

  /* start parallel region */
#pragma omp parallel default(none) private(i) \
        shared(a,b,n,stime,ftime,sum,alloctime,inittime,runtime)
  {

    /* allocate the vectors */
    #pragma omp single
    {
      stime = get_time();
      a = malloc(n*sizeof(double));
      b = malloc(n*sizeof(double));
      ftime = get_time();
      alloctime = ftime - stime;
    }

    /* initialize the vector values */
    #pragma omp master
    stime = get_time();
    #pragma omp for
    for (i=0; i<n; i++) {
      a[i] = (0.001 * (i + 1.0)) / n;
      b[i] = (0.001 * (n - i - 1.0)) / n;
    }
    #pragma omp master
    {
      ftime = get_time();
      inittime = ftime - stime;
    }

    /* start timer */
    #pragma omp master
    stime = get_time();

    /* compute dot-product */
    #pragma omp single
    sum = 0.0;
    #pragma omp for reduction(+:sum)
    for (i=0; i<n; i++)  sum += a[i]*b[i];

    /* stop timer */
    #pragma omp master
    {
      ftime = get_time();
      runtime = ftime - stime;
    }

    /* output computed value and runtime */
    #pragma omp master
    {
      printf(" vector length = %i\n",n);
#ifdef _OPENMP
      printf("   num threads = %i\n", omp_get_num_threads());
#endif
      printf("   dot-product = %.16e\n",sum);
      printf("    alloc time = %.2e\n",alloctime);
      printf("     init time = %.2e\n",inittime);
      printf("      run time = %.2e\n",runtime);
    }    

    /* free vectors */
    #pragma omp single
    {
      free(a);
      free(b);
    }

  } /* end parallel region*/

  return 0;
} /* end main */

