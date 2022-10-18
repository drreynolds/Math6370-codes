/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

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
  if (ierr != 0) {
    std::cerr << " error in MPI_Init = " << ierr << std::endl;
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != 0) {
    std::cerr << " error in MPI_Comm_size = " << ierr << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != 0) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

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

  // root node outputs parallelism information to screen
  if (myid == 0)
    std::cout << " starting MPI with " << numprocs << " processes\n";

  // determine this processor's interval
  is = ((int) (1.0*n/numprocs))*myid;
  ie = ((int) (1.0*n/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = n;

  // allocate the local portion of vectors
  stime = MPI_Wtime();
  a = new double[ie-is];
  b = new double[ie-is];
  ftime = MPI_Wtime();
  alloctime = ftime - stime;

  // initialize the vector values
  stime = MPI_Wtime();
  for (i=is; i<ie; i++) {
    a[i-is] = (0.001 * (i + 1.0)) / n;
    b[i-is] = (0.001 * (n - i - 1.0)) / n;
  }
  ftime = MPI_Wtime();
  inittime = ftime - stime;

  // start computation timer
  stime = MPI_Wtime();

  // compute dot-product
  mysum = 0.0;
  for (i=0; i<ie-is; i++)  mysum += a[i]*b[i];

  // root node collects result
  ierr = MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << " error in MPI_Reduce = " << ierr << std::endl;
    delete[] a;
    delete[] b;
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

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
