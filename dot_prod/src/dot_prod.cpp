/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   25 October 2021 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <RAJA/RAJA.hpp>

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

  // set the RAJA policies
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  std::cout << "Running dot_prod with length " << n << " with RAJA using CUDA backend:\n";

  // allocate the vectors
  std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();
  double *a, *b;
  cudaMalloc((void**)&a, n*sizeof(double));
  cudaMalloc((void**)&b, n*sizeof(double));
  std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> alloctime = ftime - stime;

  // initialize the vector values
  stime = std::chrono::system_clock::now();
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    a[i] = (0.001 * (i + 1.0)) / n;
    b[i] = (0.001 * (n - i - 1.0)) / n;
  });
  ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> inittime = ftime - stime;

  // compute dot-product
  stime = std::chrono::system_clock::now();
  RAJA::ReduceSum<rpolicy, double> sum(0.0);
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    sum += a[i]*b[i];});
  ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = ftime - stime;

  // output computed value and runtime
  std::cout << "   dot-product = " << std::setprecision(16) << sum << std::endl;
  std::cout << "    alloc time = " << alloctime.count() << std::endl;
  std::cout << "     init time = " << inittime.count() << std::endl;
  std::cout << "      run time = " << runtime.count() << std::endl;

  // delete vectors
  cudaFree(a);
  cudaFree(b);

  return 0;
} // end main
