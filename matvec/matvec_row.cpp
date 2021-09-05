/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   14 January 2013 */

// Inclusions 
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "get_time.h"


// Example routine to compute the product of an m*n matrix and an n-vector. 
int main(int argc, char* argv[]) {

  // local variables
  int m, n, i, j;
  double **A, *x, *b, runtime, norm2;
  double stime, ftime;
  std::ofstream fptr;

  // input the size of the system 
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
  A = new double*[m];
  for (i=0; i<m; i++)  A[i] = new double[n];
  x = new double[n];
  b = new double[m];
  
  // initialize the matrix and vectors 
  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      A[i][j] = 1.0/(1.0 + (i-j)*(i-j));
  for (i=0; i<m; i++)  b[i] = 0.0;
  for (j=0; j<n; j++)  x[j] = 1.0;

  // start timer 
  stime = get_time();

  // compute matrix-vector product (row-based version)
  for (i=0; i<m; i++) 
    for (j=0; j<n; j++) 
      b[i] += A[i][j]*x[j];

  // stop timer 
  ftime = get_time();
  runtime = ftime - stime;

  // output 2-norm of product and runtime to screen 
  norm2 = 0.0;
  for (i=0; i<m; i++)  norm2 += b[i]*b[i];
  std::cout << "       matrix size = " << m << " x " << n << "\n";
  std::cout << " 2-norm of product = " << sqrt(norm2) << "\n";
  std::cout << "           runtime = " << runtime << "\n";

  // output product to file 
  fptr.open("b.txt", std::fstream::out);
  for (i=0; i<m; i++)  fptr << b[i] << "\n";
  fptr.close();
  
  // free matrix and vectors 
  for (i=0; i<m; i++)  delete[] A[i];
  delete[] A;
  delete[] x;
  delete[] b;

} // end main 

