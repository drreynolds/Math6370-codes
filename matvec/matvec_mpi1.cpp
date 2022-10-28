/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "mpi.h"


// Example routine to compute the product of an m*n matrix and an n-vector.
int main(int argc, char* argv[]) {

  // intialize MPI
  int ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Init = " << ierr << "\n";
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_size = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in MPI_Comm_rank = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // root node inputs the sie of the system
  int m, n;
  if (myid == 0) {
    std::cout << "We will multiply a m*n matrix by an n-vector\n";
    std::cout << "   enter m\n";
    std::cin >> m;
    std::cout << "   enter n\n";
    std::cin >> n;
  }

  // root node sends m and n out to other processors
  ierr = MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << " error in MPI_Bcast = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << " error in MPI_Bcast = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // check for legal matrix size
  if ((m < 1) || (n < 1)) {
    if (myid == 0)
      std::cerr << " Illegal input, m = " << m << " and n = "
		<< n << " must both be >= 1\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // root node outputs parallelism information to screen
  if (myid == 0)
    std::cout << " starting MPI with " << numprocs << " processes\n";

  // allocate the matrix and vectors
  // (store entire system on every proc -- wasteful)
  int i, j;
  double **A = new double*[m];
  for (i=0; i<m; i++)  A[i] = new double[n];
  double *x = new double[n];
  double *b = new double[m];
  double *myb = new double[m];

  // initialize the matrix and vectors
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      A[i][j] = 1.0/(1.0 + (i-j)*(i-j));
  for (j=0; j<m; j++)  b[j] = 0.0;
  for (j=0; j<n; j++)  x[j] = 1.0;

  // start timer
  double stime = MPI_Wtime();

  // determine this processor's interval
  int is = ((int) (1.0*m/numprocs))*myid;
  int ie = ((int) (1.0*m/numprocs))*(myid+1);
  if (myid == numprocs-1)  ie = m;

  // compute matrix-vector product
  for (i=0; i<m; i++)  myb[i] = 0.0;
  for (i=is; i<ie; i++)
    for (j=0; j<n; j++)
      myb[i] += A[i][j]*x[j];

  // root node collects result with a reduction
  ierr = MPI_Reduce(myb, b, m, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << " error in MPI_Reduce = " << ierr << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // stop timer
  double ftime = MPI_Wtime();

  // output 2-norm of product and runtime to screen
  if (myid == 0) {
    double norm2 = 0.0;
    for (i=0; i<m; i++)  norm2 += b[i]*b[i];
    std::cout << "       matrix size = " << m << " x " << n << "\n";
    std::cout << " 2-norm of product = " << sqrt(norm2) << "\n";
    std::cout << "           runtime = " << ftime-stime << "\n";

    // output product to file
    std::ofstream fptr;
    fptr.open("b_mpi1.txt", std::fstream::out);
    for (i=0; i<m; i++)  fptr << b[i] << "\n";
    fptr.close();
  }

  // free matrix and vectors
  for (i=0; i<m; i++)  delete[] A[i];
  delete[] A;
  delete[] x;
  delete[] b;
  delete[] myb;

  // finalize MPI
  ierr = MPI_Finalize();

} // end main