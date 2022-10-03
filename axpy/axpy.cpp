/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Kokkos_Core.hpp>


// Example routine to perform some simple vector linear combinations.
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
  std::cout << "Running axpy with length " << n << " with Kokkos using OpenMP backend:\n";
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda       ExecSpace;
  typedef Kokkos::CudaSpace  MemSpace;
  std::cout << "Running axpy with length " << n << " with Kokkos using CUDA backend:\n";
#else
  typedef Kokkos::Serial     ExecSpace;
  typedef Kokkos::HostSpace  MemSpace;
  std::cout << "Running axpy with length " << n << " with Kokkos using Serial backend:\n";
#endif
  typedef Kokkos::RangePolicy<ExecSpace>    RangePol;
  typedef Kokkos::View<double*, MemSpace>   VecView;

  // Allocate the vectors as Kokkos "views"
  Kokkos::Timer timer;
  VecView x( "x", n );
  VecView y( "y", n );
  VecView z( "z", n );
  double alloctime = timer.seconds();

  // initialize a, x and y
  timer.reset();
  const double a = -3.0;
  Kokkos::parallel_for( "init_data", RangePol(0,n), KOKKOS_LAMBDA (int i) {
    x(i) = exp(2.0*(i+1)/n);
    y(i) = 1.0*(n-1)/n;
  });
  // wait for asynchronous initialization to complete, and stop initialization timer
  Kokkos::fence();
  double inittime = timer.seconds();

  // perform linear combinations
  timer.reset();
  Kokkos::parallel_for( "update_z", RangePol(0,n), KOKKOS_LAMBDA (int i) {
    z(i) = a*x(i) + y(i);
  });
  Kokkos::parallel_for( "update_x", RangePol(0,n), KOKKOS_LAMBDA (int i) {
    x(i) = y(i)/a - z(i);
  });
  Kokkos::parallel_for( "update_y", RangePol(0,n), KOKKOS_LAMBDA (int i) {
    y(i) = x(i)*y(i)/n;
  });

  // wait for asynchronous calculations to complete
  Kokkos::fence();

  // compute/output maximum value in z
  double zmax = -1e100;
  Kokkos::parallel_reduce( "max_z", RangePol(0,n), KOKKOS_LAMBDA (int i, double &mymax) {
    mymax = (mymax > z(i)) ? mymax : z(i);
  }, Kokkos::Max<double>(zmax));

  // stop timer
  double runtime = timer.seconds();

  // output computed value and runtimes
  std::cout << "      max(z) = " << std::setprecision(16) << zmax << std::endl;
  std::cout << "  alloc time = " << alloctime << std::endl;
  std::cout << "   init time = " << inittime << std::endl;
  std::cout << "    run time = " << runtime << std::endl;

  }
  Kokkos::finalize();   // note that the views x, y, and z are automatically deleted here

  return 0;
} // end main
