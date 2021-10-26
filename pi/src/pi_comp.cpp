/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   25 October 2021 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <RAJA/RAJA.hpp>


// Prototypes
RAJA_DEVICE inline double f(double a) { return (4.0 / (1.0 + a*a)); }


/* Example routine to compute pi using numerical integration via
      pi = 4 * int_0^1 1/(1+x^2) dx
   We use a simple midpoint rule for integration, over
   subintervals of fixed size 1/n, where n is a user-input parameter. */
int main(int argc, char* argv[]) {

  // ensure that an argument was passed in
  if (argc < 2) {
    std::cerr << "Error: function requires one argument (number of intervals)\n";
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  int n = atoi(argv[1]);
  if (n < 1) {
    std::cerr << "Error: number of intervals " << n << " must be greater than 0\n";
    return 1;
  }

  // set the RAJA policies
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  std::cout << "Running pi_comp with " << n << " intervals with RAJA using CUDA backend:\n";

  // start timer
  const std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();

  // set subinterval width
  const double h = 1.0 / n;

  // perform integration over n intervals
  RAJA::ReduceSum<rpolicy, double> pi(0.0);
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    double x = h * (i + 0.5);
    pi += h * f(x);});

  // stop timer
  const std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = ftime - stime;

  // output computed value and error
  std::cout << " computed pi = " << std::setprecision(16) << pi.get() << std::endl;
  std::cout << "     true pi = " << M_PI << std::endl;
  std::cout << "       error = " << M_PI-pi.get() << std::endl;
  std::cout << "     runtime = " << runtime.count() << std::endl;

  return 0;
} // end main
