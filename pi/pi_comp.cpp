/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Kokkos_Core.hpp>


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

  // initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

  // set the Kokkos execution space and range policy
#if defined(USE_OPENMP)
  typedef Kokkos::OpenMP  ExecSpace;
  std::cout << "Running pi_comp with " << n << " intervals with Kokkos using OpenMP backend:\n";
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda    ExecSpace;
  std::cout << "Running pi_comp with " << n << " intervals with Kokkos using CUDA backend:\n";
#else
  typedef Kokkos::Serial  ExecSpace;
  std::cout << "Running pi_comp with " << n << " intervals with Kokkos using Serial backend:\n";
#endif
  typedef Kokkos::RangePolicy<ExecSpace>  RangePol;

  // create/start timer
  Kokkos::Timer timer;

  // set subinterval width
  const double h = 1.0 / n;

  // perform integration over n intervals
  double pi = 0.0;
  Kokkos::parallel_reduce( RangePol(0,n), KOKKOS_LAMBDA ( int i, double &my_pi ) {
    double x = h * (i + 0.5);
    my_pi += h * (4.0 / (1.0 + x*x));
  }, pi);

  // stop timer
  double runtime = timer.seconds();

  // output computed value and error
  std::cout << " computed pi = " << std::setprecision(16) << pi << std::endl;
  std::cout << "     true pi = " << M_PI << std::endl;
  std::cout << "       error = " << M_PI-pi << std::endl;
  std::cout << "     runtime = " << runtime << std::endl;

  }
  Kokkos::finalize();

  return 0;
} // end main
