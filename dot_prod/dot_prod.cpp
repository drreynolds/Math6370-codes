/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Kokkos_Core.hpp>


// Example routine to compute the dot-product of two vectors
int main(int argc, char* argv[]) {

  // ensure that an argument was passed in
  if (argc < 2) {
    std::cerr << "Error: function requires one argument (vector length)\n";
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  int n = atoi(argv[1]);
  if (n < 1) {
    std::cerr << "Error: vector length " << n << " must be greater than 0\n";
    return 1;
  }

  // initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

  // set the Kokkos execution space, memory space, range policy, and vector view
#if defined(USE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace;
  typedef Kokkos::HostSpace  MemSpace;
  std::cout << "Running dot_prod with length " << n << " with Kokkos using OpenMP backend:\n";
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda       ExecSpace;
  typedef Kokkos::CudaSpace  MemSpace;
  std::cout << "Running dot_prod with length " << n << " with Kokkos using CUDA backend:\n";
#else
  typedef Kokkos::Serial     ExecSpace;
  typedef Kokkos::HostSpace  MemSpace;
  std::cout << "Running dot_prod with length " << n << " with Kokkos using Serial backend:\n";
#endif
  typedef Kokkos::RangePolicy<ExecSpace>    RangePol;
  typedef Kokkos::View<double*, MemSpace>   VecView;

  // Allocate the vectors as Kokkos "views"
  Kokkos::Timer timer;
  VecView a( "a", n );
  VecView b( "b", n );
  double alloctime = timer.seconds();

  // initialize the vector values on device
  timer.reset();
  Kokkos::parallel_for( RangePol(0,n), KOKKOS_LAMBDA (int i) {
    a(i) = (0.001 * (i + 1.0)) / n;
    b(i) = (0.001 * (n - i - 1.0)) / n;
  });
  // wait for asynchronous initialization to complete, and stop initialization timer
  Kokkos::fence();
  double inittime = timer.seconds();

  // compute dot-product
  timer.reset();
  double sum = 0.0;
  Kokkos::parallel_reduce( RangePol(0,n), KOKKOS_LAMBDA (int i, double &mysum) {
    mysum += a(i) * b(i);
  }, sum);
  double runtime = timer.seconds();

  // output computed value and runtimes
  std::cout << "   dot-product = " << std::setprecision(16) << sum << std::endl;
  std::cout << "    alloc time = " << alloctime << std::endl;
  std::cout << "     init time = " << inittime << std::endl;
  std::cout << "      run time = " << runtime << std::endl;

  }
  Kokkos::finalize();   // note that the views a and b are automatically deleted here

  return 0;
} // end main
