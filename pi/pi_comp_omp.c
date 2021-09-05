/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "get_time.h"


/* Prototypes */
inline double f(double a) { return (4.0 / (1.0 + a*a)); }


/* Example routine to compute pi using numerical integration via
      pi = 4 * int_0^1 1/(1+x^2) dx
   We use a simple midpoint rule for integration, over
   subintervals of fixed size 1/n, where n is a user-input parameter. */
int main(int argc, char* argv[]) {

  /* declarations */
  int i, n;
  double h, x, pi=0.0, runtime, pi_true=3.14159265358979323846;
  double stime, ftime;

  /* input the number of intervals */
  printf("Enter the number of intervals (0 quits):\n");
  i = scanf("%i", &n);
  if (n < 1 || i != 1) {
    return(-1);
  }

  /* start timer */
  stime = get_time();

  /* start OpenMP parallelism */
  #pragma omp parallel default(shared), private(x)
  {

    /* output parallelism information */
    #pragma omp single
    {
#ifdef _OPENMP
      printf(" starting OpenMP with %i processes\n", omp_get_num_threads());
#endif

      /* set subinterval width */
      h = 1.0 / n;

    } /* end omp single */

    /* perform integration over n intervals */ 
    #pragma omp for reduction(+:pi)
    for (i=0; i<n; i++) {
      x = h * (i + 0.5);
      pi += h * f(x);
    }

  } /* end omp parallel */

  /* stop timer */
  ftime = get_time();
  runtime = ftime - stime;

  /* output computed value and error */
  printf(" computed pi = %.16e\n",pi);
  printf("     true pi = %.16e\n",pi_true);
  printf("       error = %.16e\n",pi_true-pi);
  printf("     runtime = %.16e\n",runtime);

} /* end main */

