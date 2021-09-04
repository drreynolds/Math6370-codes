/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   20 January 2011 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "get_time.h"

/* Prototypes */
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // declarations
  int maxit, n, i, its;
  double lam, eps, *T, *u, *v, *w, res, runtime;
  double stime, ftime;

  // set solver input parameters
  maxit = 1000000;
  lam = 1.e-2;
  eps = 1.e-10;

  /* input the number of intervals */
  printf("Enter the number of intervals (0 quits):\n");
  i = scanf("%i", &n);
  if (n < 1 || i < 1) {
    return(-1);
  }

  /* allocate temperature and solution arrays */
  T = malloc( n * sizeof(double) );
  u = malloc( n * sizeof(double) );
  v = malloc( n * sizeof(double) );
  w = malloc( n * sizeof(double) );

  /* set random temperature field, initial guesses at chemical densities */
  for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
  for (i=0; i<n; i++)  u[i] = 0.35;
  for (i=0; i<n; i++)  v[i] = 0.1;
  for (i=0; i<n; i++)  w[i] = 0.5;

  /* start timer */
  stime = get_time();

  /* call solver over n intervals */
  for (i=0; i<n; i++) {
    chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, maxit, &its, &res);
    if (res < eps) 
      printf("    i = %i,  its = %i\n", i, its);
    else {
      printf("    error: i = %i,  its = %i,  res = %.2e,  u = %.2e,  v = %.2e,  w = %.2e\n", i, its, res, u[i], v[i], w[i]);
      return 1;
    }
  }

  /* stop timer */
  ftime = get_time();
  runtime = ftime - stime;

  /* output solution time */
  printf("     runtime = %.16e\n",runtime);

  /* free temperature and solution arrays */
  free(T);
  free(u);
  free(v);
  free(w);

  return 0;
} /* end main */

