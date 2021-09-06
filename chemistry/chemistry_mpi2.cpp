/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   4 January 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

// Prototypes
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // declarations
  int maxit, n, i, its, iend, numsent, numreceived, mycount;
  double lam, eps, *buffer, res, runtime;
  double stime, ftime;
  int ierr, numprocs, myid, sender, ansentry;
  MPI_Status status;
  bool more_work;

  // intialize MPI
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_size = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_rank = %i\n",ierr);
    return 1;
  }

  // check whether we have enough processes
  if (numprocs < 2) {
    printf("manager/worker example, requires > 1 proc! Exiting\n");
    ierr = MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }

  // set solver parameters, temporary variables
  maxit = 1000000;
  lam = 1.e-2;
  eps = 1.e-10;
  buffer = new double[3];


  /////////////////
  // manager role
  if (myid == 0) {

    printf("Starting MPI with %i processes\n\n",numprocs);

    // input the total number of intervals
    printf("Enter the total number of intervals (0 quits):\n");
    i = scanf("%i", &n);
    if (i < 1) {
      ierr = MPI_Finalize();
      return 1;
    }

    // stop for illegal n
    if (n < 1) {
      ierr = MPI_Finalize();
      return 1;
    }

    // allocate temperature field and solution arrays
    double *T = new double[n];
    double *u = new double[n];
    double *v = new double[n];
    double *w = new double[n];

    // set random temperature field
    for (i=0; i<n; i++)  
      T[i] = random() / (pow(2.0,31.0) - 1.0);

    // start timer
    stime = MPI_Wtime();

    // do manager/worker computation (assumes n >= numprocs)
    numsent = 0;

    // first send a part to each worker
    iend = (n < numprocs-1) ? n : numprocs-1;
    for (i=0; i<iend; i++) {
      // fill send buffer
      buffer[0] = T[i];
      // send with tag as entry in temperature array
      ierr = MPI_Send(buffer, 1, MPI_DOUBLE, i+1, numsent, MPI_COMM_WORLD);
      if (ierr != 0) {
	printf("Error in MPI_Send = %i\n",ierr);
	ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	return 1;
      }
      numsent++;
    } // for i

    // receive work from each proc, and assign additional work if available
    for (i=0; i<n; i++) {

      // receive answers from any process
      ierr = MPI_Recv(buffer, 3, MPI_DOUBLE, MPI_ANY_SOURCE, 
		      MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (ierr != 0) {
	printf("Error in MPI_Recv = %i\n",ierr);
	ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	return 1;
      }
      // decode the sender and solution entry from status
      sender   = status.MPI_SOURCE;
      ansentry = status.MPI_TAG;
      // store results
      u[ansentry] = buffer[0];
      v[ansentry] = buffer[1];
      w[ansentry] = buffer[2];
      if (numsent < n) {  // send another row
	buffer[0] = T[numsent];
	ierr = MPI_Send(buffer, 1, MPI_DOUBLE, sender, numsent, MPI_COMM_WORLD);
	if (ierr != 0) {
	  printf("Error in MPI_Send = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	  return 1;
	}
	numsent++;
      } else {  // tell senders that work is complete, by sending message of zero size
	ierr = MPI_Send(buffer, 0, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD);
	if (ierr != 0) {
	  printf("Error in MPI_Send = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	  return 1;
	}
      } // if numsent
           
    } // for i

    // stop timer
    ftime = MPI_Wtime();
    runtime = ftime - stime;

    // output runtime
    printf("     runtime = %.16e\n",runtime);

    // delete temperature and solution arrays
    delete[] T;
    delete[] u;
    delete[] v;
    delete[] w;


  /////////////////
  // worker role
  } else {
    
    // allocate temperature and solution values to reuse
    double T, u, v, w;

    // initialize internal work counter
    mycount = 0;

    // loop over assignments from manager
    more_work = true;
    while (more_work) {

      // receive message from manager
      ierr = MPI_Recv(buffer, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, 
		      MPI_COMM_WORLD, &status);
      if (ierr != 0) {
	printf("Error in MPI_Recv = %i\n",ierr);
	ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	return 1;
      }

      // check if work is complete
      ierr = MPI_Get_count(&status, MPI_DOUBLE, &numreceived);
      if (numreceived == 0) {
	more_work = false;

      // otherwise, do computation
      } else {

	// set row index
	i = status.MPI_TAG;

	// set input and initial guess
	T = buffer[0];
	u = 0.35;
	v = 0.1;
	w = 0.5;

	// call solver
	chem_solver(T, &u, &v, &w, lam, eps, maxit, &its, &res);
	if (res < eps) {
 	  // printf("    i = %i,  its = %i\n", i, its);
	}
	else {
	  printf("    error: i=%i, its=%i, res=%.2e, u=%.2e, v=%.2e, w=%.2e\n", 
		 i, its, res, u, v, w);
	  ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	  return 1;
	}

	// send result to manager
	buffer[0] = u;
	buffer[1] = v;
	buffer[2] = w;
	ierr = MPI_Send(buffer, 3, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
	if (ierr != 0) {
	  printf("Error in MPI_Send = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, 1);
	  return 1;
	}

	// update work counter
	mycount++;

      } // if status.MPI_TAG

    } // end while loop

    // output work performed
    printf(" proc %i computed %i iterations\n", myid, mycount);

  } // end worker section

  // finalize MPI
  ierr = MPI_Finalize();

} // end main
