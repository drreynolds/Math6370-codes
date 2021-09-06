/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "mpi.h"


// Example routine to compute the dot-product of two vectors
int main(int argc, char* argv[]) {

  // declarations
  int i, n, is, ie;
  double *a, *b, sum, mysum, alloctime, inittime, runtime;
  double stime, ftime;
  int ierr, numprocs, myid;

  // initialize MPI
  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  // ensure that an argument was passed in
  if (argc < 2) {
    std::cerr << "Error: function requires one argument (vector length)\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // set n as the input argument, and ensure it's positive
  n = atoi(argv[1]);
  if (n < 1) {
    std::cerr << "Error: vector length " << n << " must be greater than 0\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // allocate the vectors (go ahead and put the whole thing on all procs)
  stime = MPI_Wtime();
  a = new double[n];
  b = new double[n];
  ftime = MPI_Wtime();
  alloctime = ftime - stime;

  // initialize the vector values
  stime = MPI_Wtime();
  for (i=0; i<n; i++) {
    a[i] = (0.001 * (i + 1.0)) / n;
    b[i] = (0.001 * (n - i - 1.0)) / n;
  }
  ftime = MPI_Wtime();
  inittime = ftime - stime;

  // start computation timer
  stime = MPI_Wtime();

  // determine this processor's interval
  is = ((int) (1.0*n/numprocs))*myid;
  ie = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = n;

  // compute dot-product
  mysum = 0.0;
  for (i=is; i<ie; i++)  mysum += a[i]*b[i];

  // root node collects result
  ierr = MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // stop timer
  ftime = MPI_Wtime();
  runtime = ftime-stime;

  // output computed value and runtime
  if (myid == 0) {
    std::cout << " vector length = " << n << std::endl;
    std::cout << "   dot-product = " << std::setprecision(16) << sum << std::endl;
    std::cout << "    alloc time = " << std::setprecision(2) << alloctime << std::endl;
    std::cout << "     init time = " << inittime << std::endl;
    std::cout << "      run time = " << runtime << std::endl;
  }

  // delete vectors
  delete[] a;
  delete[] b;

  // finalize MPI
  ierr = MPI_Finalize();

} // end main

