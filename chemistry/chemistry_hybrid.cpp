/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   19 April 2017 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "mpi.h"
#ifdef _OPENMP
#include "omp.h"
#endif


// Prototypes
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // intialize MPI (we don't actually need 'MPI_THREAD_MULTIPLE', I just want 
  // to inquire about the maximum level supported in this MPI implementation
  int ierr, numprocs, myid, provided;
  ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Init = " << ierr << std::endl;
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_size = " << ierr << std::endl;
    return 1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    std::cerr << "Error in MPI_Comm_rank = " << ierr << std::endl;
    return 1;
  }

  // set solver input parameters
  int maxit = 1000000;
  double lam = 1.e-2;
  double eps = 1.e-10;

  // root node outputs parallelism information to screen
  if (myid == 0) {
    std::cout << "Starting MPI with " << numprocs << " processes\n";
    switch(provided) {
    case MPI_THREAD_SINGLE:
      std::cout << "MPI library supports threading at 'single' level\n";
      break;
    case MPI_THREAD_FUNNELED:
      std::cout << "MPI library supports threading at 'funneled' level\n";
      break;
    case MPI_THREAD_SERIALIZED:
      std::cout << "MPI library supports threading at 'serialized' level\n";
      break;
    case MPI_THREAD_MULTIPLE:
      std::cout << "MPI library supports threading at 'multiple' level\n";
      break;
    }
  }
  
  // ensure that an argument was passed in
  if (argc < 2) {
    if (myid == 0)
      std::cerr << "Error: function requires one argument (number of intervals)\n";
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  int n = atoi(argv[1]);
  if (n < 1) {
    if (myid == 0)
      std::cerr << "Error: number of intervals = " << n << " (must be greater than 0)\n";
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // determine this processor's interval
  int is = ((int) 1.0*n/numprocs)*myid;
  int ie = ((int) 1.0*n/numprocs)*(myid+1) - 1;
  if (myid == numprocs-1)  ie = n-1;

  // declare data arrays
  double *allT=NULL, *allu=NULL, *allv=NULL, *allw=NULL, *T=NULL, *u=NULL, *v=NULL, *w=NULL;

  // root node allocates full temperature field and solution arrays
  if (myid == 0) {
    allT = new double[n];
    allu = new double[n];
    allv = new double[n];
    allw = new double[n];

    // set random temperature field
    for (int i=0; i<n; i++)  allT[i] = random() / (pow(2.0,31.0) - 1.0);
  }
  
  // all procs allocate local temperature field and solution arrays
  T = new double[ie-is+1];
  u = new double[ie-is+1];
  v = new double[ie-is+1];
  w = new double[ie-is+1];

  // set initial guesses at chemical densities
  for (int i=0; i<(ie-is+1); i++)  u[i] = 0.35;
  for (int i=0; i<(ie-is+1); i++)  v[i] = 0.1;
  for (int i=0; i<(ie-is+1); i++)  w[i] = 0.5;

  // allocate gatherv/scatterv temporary arrays
  int *counts = new int[numprocs];
  int *displs = new int[numprocs];

  // start timer
  double stime = MPI_Wtime();

  // root node sends out portions of temperature field to all procs
  for (int i=0; i<numprocs; i++) {
    int js = floor(1.0*n/numprocs)*i;
    int je = floor(1.0*n/numprocs)*(i+1) - 1;
    if (i == numprocs-1)  je = n-1;
    counts[i] = je-js+1;
    displs[i] = js;
  }
  ierr = MPI_Scatterv(allT, counts, displs, MPI_DOUBLE, T, 
		      ie-is+1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  double commtime = MPI_Wtime()-stime;
  
  // declare shared variables for inside OpenMP region
  int proccount=0, err=0;

  // start OpenMP parallel region for actual calculations
  #pragma omp parallel default(shared)
  {

    // declare private data for thread
    int numthreads=1;
    int tid=0;

    // master thread on each MPI task outputs threading information
#ifdef _OPENMP
    numthreads = omp_get_num_threads();
    tid = omp_get_thread_num();
    #pragma omp master
    std::cout << "MPI task " << myid << " starting OpenMP with " << numthreads << " threads\n";
#endif

    // start thread timer
    double threadstart = MPI_Wtime();
    
    // MPI tasks call solver over local intervals
    int threadcount=0;
    #pragma omp for schedule(guided) reduction(+:err)
    for (int i=0; i<(ie-is+1); i++) {
      int its;
      double res;
      chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, maxit, &its, &res);
      if (res > eps) {
        std::cout << "    error: proc=" << myid << ", thread=" << tid << ", i=" << is+i << ", its=" << its
                  << ", res=" << res << ", u=" << u[i] << ", v=" << v[i] << ", w=" << w[i] << std::endl;
        err += 1;
      }

      // update thread increment
      threadcount++;
      
    } // end for i (end omp for)

    // threads output loop partitioning information
    #pragma omp critical
    {
      double threadtime = MPI_Wtime() - threadstart;
      std::cout << " proc " << myid << " / thread " << tid << " computed " << threadcount << " iterations in "
                << threadtime << " sec.\n";
      proccount += threadcount;
    }

  } // end omp parallel

  if (err > 0) {
    std::cout << "*** PROC " << myid << " WARNING ***\n" << "  See above for solver error messages.\n";
  }
  
  // root node collects results
  double stime2 = MPI_Wtime();
  ierr = MPI_Gatherv(u, ie-is+1, MPI_DOUBLE, allu, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << "Error in MPI_Gatherv = " << ierr << std::endl;
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  ierr = MPI_Gatherv(v, ie-is+1, MPI_DOUBLE, allv, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << "Error in MPI_Gatherv = " << ierr << std::endl;
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  ierr = MPI_Gatherv(w, ie-is+1, MPI_DOUBLE, allw, counts, 
		     displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (ierr != 0) {
    std::cerr << "Error in MPI_Gatherv = " << ierr << std::endl;
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // stop timer
  double ftime = MPI_Wtime();
  double runtime = ftime - stime;
  commtime += ftime-stime2;

  // output solution time on each proc
  std::cout << " proc " << myid << " computed " << proccount << " total iterations, runtime = " << runtime
            << ", comm time = " << commtime << std::endl;

  // free arrays
  delete[] T;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] counts;
  delete[] displs;
  if (allT != NULL) delete[] allT;
  if (allu != NULL) delete[] allu;
  if (allv != NULL) delete[] allv;
  if (allw != NULL) delete[] allw;

  
  // finalize MPI
  ierr = MPI_Finalize();

} // end main
