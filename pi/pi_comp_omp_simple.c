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

  /* local variables */
  int n, i;
  double h, x, pi, runtime, pi_true=3.14159265358979323846;
  double stime, ftime;

  /* input the number of intervals */
  printf("Enter the number of intervals (0 quits):\n");
  i = scanf("%i", &n);
  if (n < 1 || i != 1) {
    return(-1);
  }

  /* start timer */
  stime = get_time();

  /* set subinterval width */
  h = 1.0 / n;

  /* perform integration over n intervals */ 
  pi = 0.0;
#pragma omp parallel for reduction(+:pi) private(x)
  for (i=0; i<n; i++) {
    x = h * (i + 0.5);
    pi += h * f(x);
  }

  /* stop timer */
  ftime = get_time();
  runtime = ftime - stime;

  /* output computed value and error */
  printf(" computed pi = %.16e\n",pi);
  printf("     true pi = %.16e\n",pi_true);
  printf("       error = %.16e\n",pi_true-pi);
  printf("     runtime = %.16e\n",runtime);

} /* end main */
