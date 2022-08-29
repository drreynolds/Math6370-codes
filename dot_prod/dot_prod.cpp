/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include "get_time.h"


// Example routine to compute the dot-product of two vectors
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
  double stime = get_time();
  double *a = new double[n];
  double *b = new double[n];
  double ftime = get_time();
  double alloctime = ftime-stime;

  // initialize the vector values
  stime = get_time();
  for (int i=0; i<n; i++)    a[i] = (0.001 * (i + 1.0)) / n;
  for (int i=0; i<n; i++)    b[i] = (0.001 * (n - i - 1.0)) / n;
  ftime = get_time();
  double inittime = ftime-stime;

  // compute dot-product
  stime = get_time();
  double sum = 0.0;
  for (int i=0; i<n; i++)   sum += a[i]*b[i];
  ftime = get_time();
  double runtime = ftime - stime;

  // output computed value and runtime
  printf(" vector length = %i\n",n);
  printf("   dot-product = %.16e\n",sum);
  printf("    alloc time = %.2e\n",alloctime);
  printf("     init time = %.2e\n",inittime);
  printf("      run time = %.2e\n",runtime);

  // delete vectors
  delete[] a;
  delete[] b;

  return 0;
} // end main
