/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   13 January 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_time.h"
#include <RAJA/RAJA.hpp>


// Example routine to perform some simple vector linear combinations.
int main(int argc, char* argv[]) {

  // ensure that an argument was passed in
  if (argc < 2) {
    printf("Error: function requires one argument (vector length)\n");
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  int n = atoi(argv[1]);
  if (n < 1) {
    printf("Error: vector length %i must be greater than 0\n", n);
    return 1;
  }

  // allocate the vectors
  double *x, *y, *z;
#if defined(RAJA_ENABLE_CUDA)
  cudaMallocManaged((void**)&x, n*sizeof(double), cudaMemAttachGlobal);
  cudaMallocManaged((void**)&y, n*sizeof(double), cudaMemAttachGlobal);
  cudaMallocManaged((void**)&z, n*sizeof(double), cudaMemAttachGlobal);
#else
  x = new double[n];
  y = new double[n];
  z = new double[n];
#endif

  // set the RAJA policies
#if defined(RAJA_ENABLE_CUDA)
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  printf("Running axpy with length %i with RAJA using CUDA backend:\n", n);
#elif defined(RAJA_ENABLE_OPENMP)
  using epolicy = RAJA::omp_parallel_for_exec;
  using rpolicy = RAJA::omp_reduce;
  printf("Running axpy with length %i with RAJA using OpenMP backend:\n", n);
#else
  using epolicy = RAJA::seq_exec;
  using rpolicy = RAJA::seq_reduce;
  printf("Running axpy with length %i with RAJA using sequential backend:\n", n);
#endif


  // start timer
  double stime = get_time();

  // initialize a, x and y
  double a = -3.0;
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_HOST_DEVICE (int i) {
    x[i] = exp(2.0*(i+1)/n);
    y[i] = 1.0*(n-1)/n;
  });

  // perform linear combinations
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_HOST_DEVICE (int i) {
    z[i] = a*x[i] + y[i];});
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_HOST_DEVICE (int i) {
    x[i] = y[i]/a - z[i];});
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_HOST_DEVICE (int i) {
    y[i] = x[i]*y[i]/n;});

  // output maximum value in z
  RAJA::ReduceMax<rpolicy, double> zmax(-1e100);
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_HOST_DEVICE (int i) {
    zmax.max(z[i]);});
  printf("  max(z) = %.16e\n",zmax.get());

  // stop timer
  double ftime = get_time();
  double runtime = ftime - stime;

  // output total time
  printf(" runtime = %.16e\n",runtime);

  // free vectors
#if defined(RAJA_ENABLE_CUDA)
  cudaFree(x);
  cudaFree(y);
  cudaFree(z);
#else
  delete[] x;
  delete[] y;
  delete[] z;
#endif

  return 0;
} // end main
