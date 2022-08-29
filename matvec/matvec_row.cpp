/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "get_time.h"


// Example routine to compute the product of an m*n matrix and an n-vector.
int main(int argc, char* argv[]) {

  // input the size of the system
  int m, n;
  std::cout << "We will multiply a m*n matrix by an n-vector\n";
  std::cout << "   enter m\n";
  std::cin >> m;
  std::cout << "   enter n\n";
  std::cin >> n;
  if ((m < 1) || (n < 1)) {
    std::cerr << " Illegal input, m = " << m << " and n = "
	      << n << " must both be >= 1\n";
    return 1;
  }

  // allocate the matrix and vectors
  int i, j;
  double **A = new double*[m];
  for (i=0; i<m; i++)  A[i] = new double[n];
  double *x = new double[n];
  double *b = new double[m];

  // initialize the matrix and vectors
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      A[i][j] = 1.0/(1.0 + (i-j)*(i-j));
  for (i=0; i<m; i++)  b[i] = 0.0;
  for (j=0; j<n; j++)  x[j] = 1.0;

  // start timer
  double stime = get_time();

  // compute matrix-vector product (row-based version)
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      b[i] += A[i][j]*x[j];

  // stop timer
  double ftime = get_time();
  double runtime = ftime - stime;

  // output 2-norm of product and runtime to screen
  double norm2 = 0.0;
  for (i=0; i<m; i++)  norm2 += b[i]*b[i];
  std::cout << "       matrix size = " << m << " x " << n << "\n";
  std::cout << " 2-norm of product = " << sqrt(norm2) << "\n";
  std::cout << "           runtime = " << runtime << "\n";

  // output product to file
  std::ofstream fptr;
  fptr.open("b.txt", std::fstream::out);
  for (i=0; i<m; i++)  fptr << b[i] << "\n";
  fptr.close();

  // free matrix and vectors
  for (i=0; i<m; i++)  delete[] A[i];
  delete[] A;
  delete[] x;
  delete[] b;

} // end main
