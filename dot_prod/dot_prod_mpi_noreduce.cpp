/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions 
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include "mpi.h"


// Example routine to compute the dot-product of two vectors. 
int main(int argc, char* argv[]) {

  // declarations
  int i, n, p, is, ie;
  double *a, *b, sum, mysum, alloctime, inittime, runtime;
  double stime, ftime;
  int numprocs, myid;
  MPI_Status stat;

  // intialize MPI 
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Init\n";
    return 1;
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &numprocs) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_size\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if (MPI_Comm_rank(MPI_COMM_WORLD, &myid) != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_rank\n";
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

  // perform manual reduction
  if (myid != 0) {
    // everyone sends value to root proc
    if (MPI_Send(&mysum, 1, MPI_DOUBLE, 0, 
		 myid, MPI_COMM_WORLD) != MPI_SUCCESS) {
      std::cerr << "Error in MPI_Send\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  } else {
    // root receives values from others and adds to own
    sum = mysum;
    for (p=1; p<numprocs; p++) {
      if (MPI_Recv(&mysum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, 
		   MPI_COMM_WORLD, &stat) != MPI_SUCCESS) {
	std::cerr << "Error in MPI_Send\n";
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
      sum += mysum;
    }
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

  // free vectors 
  delete[] a;
  delete[] b;

  // finalize MPI 
  MPI_Finalize();

} // end main 

