/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "get_time.h"

// Prototypes
void chem_solver(double, double*, double*, double*,
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // set solver input parameters
  int maxit = 1000000;
  double lam = 1.e-2;
  double eps = 1.e-10;

  // input the number of intervals
  int n;
  printf("Enter the number of intervals (0 quits):\n");
  scanf("%i", &n);
  if (n < 1)  return 1;

  // allocate temperature and solution arrays
  double *T = new double[n];
  double *u = new double[n];
  double *v = new double[n];
  double *w = new double[n];

  // set random temperature field, initial guesses at chemical densities
  int i;
  for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
  for (i=0; i<n; i++)  u[i] = 0.35;
  for (i=0; i<n; i++)  v[i] = 0.1;
  for (i=0; i<n; i++)  w[i] = 0.5;

  // start timer
  double stime = get_time();

  // call solver over n intervals
  for (i=0; i<n; i++) {
    int its;
    double res;
    chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, maxit, &its, &res);
    if (res < eps)
      printf("    i = %i,  its = %i\n", i, its);
    else {
      printf("    error: i=%i, its=%i, res=%.2e, u=%.2e, v=%.2e, w=%.2e\n",
	     i, its, res, u[i], v[i], w[i]);
      return 1;
    }
  }

  // stop timer
  double ftime = get_time();
  double runtime = ftime - stime;

  // output solution time
  printf("     runtime = %.16e\n",runtime);

  // delete temperature and solution arrays
  delete[] T;
  delete[] u;
  delete[] v;
  delete[] w;

  return 0;
} // end main
