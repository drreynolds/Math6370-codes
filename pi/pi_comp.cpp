/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "get_time.h"


// Prototypes
inline double f(double a) { return (4.0 / (1.0 + a*a)); }


/* Example routine to compute pi using numerical integration via
      pi = 4 * int_0^1 1/(1+x^2) dx
   We use a simple midpoint rule for integration, over
   subintervals of fixed size 1/n, where n is a user-input parameter. */
int main(int argc, char* argv[]) {

  // declarations
  int i, n;
  double h, x, pi, runtime, pi_true=3.14159265358979323846;
  double stime, ftime;

  // input the number of intervals
  std::cout << "Enter the number of intervals (0 quits):\n";
  std::cin >> n;
  if (n<1) {
    return 1;
  }

  // start timer
  stime = get_time();

  // set subinterval width
  h = 1.0 / n;

  // perform integration over n intervals
  pi = 0.0;
  for (i=0; i<n; i++) {
    x = h * (i + 0.5);
    pi += h * f(x);
  }

  // stop timer
  ftime = get_time();
  runtime = ftime - stime;

  // output computed value and error
  std::cout << " computed pi = " << std::setprecision(16) << pi << std::endl;
  std::cout << "     true pi = " << pi_true << std::endl;
  std::cout << "       error = " << pi_true-pi << std::endl;
  std::cout << "     runtime = " << runtime << std::endl;

  return 0;
} // end main

