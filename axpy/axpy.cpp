/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions 
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include "get_time.h"


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
  double *x = new double[n];
  double *y = new double[n];
  double *z = new double[n];

  // start timer 
  double stime = get_time();

  // initialize a, x and y 
  double a = -3.0;
  int i;
  for (i=0; i<n; i++) {
    x[i] = exp(2.0*(i+1)/n);
    y[i] = 1.0*(n-1)/n;
  }

  // perform linear combinations 
  for (i=0; i<n; i++)  z[i] = a*x[i] + y[i];
  for (i=0; i<n; i++)  x[i] = y[i]/a - z[i];
  for (i=0; i<n; i++)  y[i] = x[i]*y[i]/n;

  // output maximum value in z 
  double zmax = z[0];
  for (i=0; i<n; i++)  zmax = (zmax > z[i]) ? zmax : z[i];
  std::cout << "  max(z) = " << std::setprecision(16) << zmax << std::endl;

  // stop timer 
  double ftime = get_time();
  double runtime = ftime - stime;

  // output total time 
  std::cout << " runtime = " << std::setprecision(16) << runtime << std::endl;

  // free vectors 
  delete[] x;
  delete[] y;
  delete[] z;

  return 0;
} // end main 

