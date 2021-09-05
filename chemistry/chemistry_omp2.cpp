/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   14 January 2013 */


// Inclusions 
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "get_time.h"


// Prototypes 
int chem_solver_orphan(int, double*, double*, double*, 
		       double*, double, double, int, int*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver_orphan.c */
int main(int argc, char* argv[]) {

  // declarations
  int maxit, n, i, its, mycount, numprocs, err=0, locerr;
  double lam, eps, *T, *u, *v, *w, res, stime, ftime;

  // set solver input parameters
  maxit = 1000000;
  lam = 1.e-2;
  eps = 1.e-10;

  // input the number of intervals
  std::cout << "Enter the number of intervals (0 quits):\n";
  std::cin >> n;
  if (n < 1) {
    return 1;
  }

  // allocate temperature and solution arrays
  T = new double[n];
  u = new double[n];
  v = new double[n];
  w = new double[n];

  // set random temperature field, initial guesses at chemical densities 
  for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);
  for (i=0; i<n; i++)  u[i] = 0.35;
  for (i=0; i<n; i++)  v[i] = 0.1;
  for (i=0; i<n; i++)  w[i] = 0.5;

  // start timer 
  stime = get_time();

  // start OpenMP parallelism 
  #pragma omp parallel default(shared) private(mycount,err)
  {

    // output parallelism information 
#ifdef _OPENMP
    #pragma omp master
    {
      numprocs = omp_get_num_threads();
      std::cout << " starting OpenMP with " << numprocs << " processes\n";
    }
#endif    

    // call solver 
    locerr = chem_solver_orphan(n, T, u, v, w, lam, eps, maxit, &mycount);

    // accumulate errors 
    #pragma omp critical
    err += locerr;

    // output loop partitioning information 
#ifdef _OPENMP
    #pragma omp critical    
    std::cout << " thread " << omp_get_thread_num() << " computed " 
	      << mycount << " iterations\n";
#endif

  } // end omp parallel 
  if (err > 0) {
    std::cout << "*** WARNING ***\n";
    std::cout << "  See above messages for errors in solver\n";
  }

  // stop timer 
  ftime = get_time();

  // output solution time 
  std::cout << "     runtime = " << ftime - stime << std::endl;

  // free temperature and solution arrays 
  delete[] T;
  delete[] u;
  delete[] v;
  delete[] w;

  return 0;
} // end main 

