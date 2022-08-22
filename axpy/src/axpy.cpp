/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   25 October 2021 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <RAJA/RAJA.hpp>


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

  // allocate the vectors
  double *x, *y, *z;
  cudaMalloc((void**)&x, n*sizeof(double));
  cudaMalloc((void**)&y, n*sizeof(double));
  cudaMalloc((void**)&z, n*sizeof(double));

  // set the RAJA policies
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  std::cout << "Running axpy with length " << n << " with RAJA using CUDA backend:\n";

  // start timer
  const std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();

  // initialize a, x and y
  const double a = -3.0;
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    x[i] = exp(2.0*(i+1)/n);
    y[i] = 1.0*(n-1)/n;
  });

  // perform linear combinations
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    z[i] = a*x[i] + y[i];});
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    x[i] = y[i]/a - z[i];});
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    y[i] = x[i]*y[i]/n;});

  // output maximum value in z
  RAJA::ReduceMax<rpolicy, double> zmax(-1e100);
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    zmax.max(z[i]);});
  std::cout << "  max(z) = " << std::setprecision(16) << zmax.get() << std::endl;

  // stop timer
  const std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = ftime - stime;

  // output total time
  std::cout << " runtime = " << std::setprecision(16) << runtime.count() << std::endl;

  // free vectors
  cudaFree(x);
  cudaFree(y);
  cudaFree(z);

  return 0;
} // end main