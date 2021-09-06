/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   20 January 2011 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "mpi.h"


/* Prototypes */
void chem_solver(double, double*, double*, double*, 
		 double, double, int, int*, double*);


/* Example routine to compute the equilibrium chemical densities at 
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the 
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  /* intialize MPI */
  int ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Init = %i\n",ierr);
    return -1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_size = %i\n",ierr);
    return -1;
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    printf("Error in MPI_Comm_rank = %i\n",ierr);
    return -1;
  }

  /* check whether we have enough processors */
  if (numprocs < 2) {
    printf("master/worker example, requires > 1 proc! Exiting\n");
    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
    return -1;
  }

  /* set solver parameters, temporary variables */
  int chunk = 10;
  int maxit = 1000000;
  double lam = 1.e-2;
  double eps = 1.e-10;
  double *buffer = malloc( 3 * chunk * sizeof(double) );
  MPI_Status status;
  int i, sender, ansentry;


  /***************/
  /* master role */
  if (myid == 0) {

    /* input the total number of intervals */
    int n;
    printf("Enter the total number of intervals (0 quits):\n");
    i = scanf("%i", &n);
    if (i < 1) {
      ierr = MPI_Finalize();
      return -1;
    }
    printf(" starting MPI with %i processes\n",numprocs);

    /* stop for illegal n */
    if (n < 1) {
      ierr = MPI_Finalize();
      return -1;
    }

    /* allocate temperature field and solution arrays */
    double *T = malloc( n * sizeof(double) );
    double *u = malloc( n * sizeof(double) );
    double *v = malloc( n * sizeof(double) );
    double *w = malloc( n * sizeof(double) );

    /* set random temperature field */
    for (i=0; i<n; i++)  T[i] = random() / (pow(2.0,31.0) - 1.0);

    /* start timer */
    time_t stime = time(NULL);

    /* do master/slave computation (assumes n >= numprocs) */
    int numsent = 0;

    /* first send a part to each slave */
    for (i=0; i<numprocs-1; i++) {
      /* fill send buffer */
      if (numsent + chunk <= n)
	for (i=0; i<chunk; i++) buffer[i] = T[numsent+i];
      else {
	for (i=0; i<n-numsent; i++) buffer[i] = T[numsent+i];
	for (i=n=numsent; i<chunk; i++) buffer[i] = 1.0;
      }
      /* send with tag as entry */
      ierr = MPI_Send(buffer, chunk, MPI_DOUBLE, i+1, i, MPI_COMM_WORLD);
      if (ierr != 0) {
	printf("Error in MPI_Send = %i\n",ierr);
	ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	return -1;
      }
      numsent += chunk;
    } /* for i */

    /* if all work was done in first phase, receive work only */
    if (numsent >= n-1) {

      /* receive work from each proc */
      for (i=1; i<numprocs; i++) {

	/* receive answers from any processor */
	ierr = MPI_Recv(buffer, 3*chunk, MPI_DOUBLE, MPI_ANY_SOURCE, 
			MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (ierr != 0) {
	  printf("Error in MPI_Recv = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	  return -1;
	}
	/* decode the sender and solution entry from status */
	sender = status.MPI_SOURCE;
	ansentry = status.MPI_TAG;
	/* store results */
	if (ansentry+chunk-1 <= n) {
	  for (i=0; i<chunk; i++) u[ansentry+i] = buffer[i];
	  for (i=0; i<chunk; i++) v[ansentry+i] = buffer[chunk+i];
	  for (i=0; i<chunk; i++) w[ansentry+i] = buffer[2*chunk+i];
	} else {
	  for (i=0; i<n-ansentry+1; i++)  u[ansentry+i] = buffer[i];
	  for (i=0; i<n-ansentry+1; i++)  v[ansentry+i] = buffer[chunk+i];
	  for (i=0; i<n-ansentry+1; i++)  w[ansentry+i] = buffer[2*chunk+i];
	}
	
	/* tell senders that work is complete */
	ierr = MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, -1, MPI_COMM_WORLD);
	if (ierr != 0) {
	  printf("Error in MPI_Send = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	  return -1;
	}
	
      } /* end for i */
      
    } else {  

      /* more work left; receive work from each proc, assign more if available*/
      for (i=1; i<=ceil(1.0*n/chunk); i++) {
	
	/* receive answers from any processor */
	ierr = MPI_Recv(buffer, 3*chunk, MPI_DOUBLE, MPI_ANY_SOURCE, 
			MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	if (ierr != 0) {
	  printf("Error in MPI_Recv = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	  return -1;
	}
	/* decode the sender and solution entry from status */
	sender = status.MPI_SOURCE;
	ansentry = status.MPI_TAG;
	/* store results */
	if (ansentry+chunk-1 <= n) {
	  for (i=0; i<chunk; i++) u[ansentry+i] = buffer[i];
	  for (i=0; i<chunk; i++) v[ansentry+i] = buffer[chunk+i];
	  for (i=0; i<chunk; i++) w[ansentry+i] = buffer[2*chunk+i];
	} else {
	  for (i=0; i<n-ansentry+1; i++)  u[ansentry+i] = buffer[i];
	  for (i=0; i<n-ansentry+1; i++)  v[ansentry+i] = buffer[chunk+i];
	  for (i=0; i<n-ansentry+1; i++)  w[ansentry+i] = buffer[2*chunk+i];
	}
	if (numsent < n) {  /* send another set of rows, if available */
	  if (numsent + chunk <= n)
	    for (i=0; i<chunk; i++) buffer[i] = T[numsent+i];
	  else {
	    for (i=0; i<n-numsent; i++) buffer[i] = T[numsent+i];
	    for (i=n=numsent; i<chunk; i++) buffer[i] = 1.0;
	  }
	  ierr = MPI_Send(buffer, chunk, MPI_DOUBLE, sender, numsent, 
			  MPI_COMM_WORLD);
	  if (ierr != 0) {
	    printf("Error in MPI_Send = %i\n",ierr);
	    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	    return -1;
	  }
	  numsent += chunk;
	} else {  /* tell senders that work is complete */
	  ierr = MPI_Send(MPI_BOTTOM, 0, MPI_DOUBLE, sender, -1, MPI_COMM_WORLD);
	  if (ierr != 0) {
	    printf("Error in MPI_Send = %i\n",ierr);
	    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	    return -1;
	  }
	} /* if numsent */
           
      } /* for i */

    } /* if numsent */

    /* stop timer */
    time_t ftime = time(NULL);
    double runtime = ((double) (ftime - stime));

    /* output runtime */
    printf("     runtime = %.16e\n",runtime);

    /* free temperature and solution arrays */
    free(T);
    free(u);
    free(v);
    free(w);

  } else {    /* end master section, start worker section */
    
    /* allocate temperature field and solution arrays */
    double *T = malloc( chunk * sizeof(double) );
    double *u = malloc( chunk * sizeof(double) );
    double *v = malloc( chunk * sizeof(double) );
    double *w = malloc( chunk * sizeof(double) );

    /* initialize internal work counter */
    int mycount=0;

    /* loop over assignments from master */
    int more_work = 1;
    int ibase;
    while (more_work) {

      /* receive message from master */
      ierr = MPI_Recv(buffer, chunk, MPI_DOUBLE, 0, MPI_ANY_TAG, 
		      MPI_COMM_WORLD, &status);
      if (ierr != 0) {
	printf("Error in MPI_Recv = %i\n",ierr);
	ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	return -1;
      }

      /* check if work is complete */
      if (status.MPI_TAG == -1)  
	more_work = 0;

      /* otherwise, do computation */
      else {

	/* set row index */
	ibase = status.MPI_TAG;

	/* set input and initial guess */
	for (i=0; i<chunk; i++)  T[i] = buffer[i];
	for (i=0; i<chunk; i++)  u[i] = 0.35;
	for (i=0; i<chunk; i++)  v[i] = 0.1;
	for (i=0; i<chunk; i++)  w[i] = 0.5;

	/* call solver over entries */
	int its;
	double res;
	for (i=0; i<chunk; i++) {
	  chem_solver(T[i], &(u[i]), &(v[i]), &(w[i]), lam, eps, 
		      maxit, &its, &res);
	  if (res < eps) {
/* 	    printf("    i = %i,  its = %i\n", ibase+i, its);  */
	  }
	  else {
	    printf("    error: i=%i, its=%i, res=%.2e, u=%.2e, v=%.2e, w=%.2e\n", 
		   ibase+i, its, res, u[i], v[i], w[i]);
	    ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	    return -1;
	  }
	}

	/* send results to master */
	for (i=0; i<chunk; i++)  buffer[i] = u[i];
	for (i=0; i<chunk; i++)  buffer[chunk+i] = v[i];
	for (i=0; i<chunk; i++)  buffer[2*chunk+i] = w[i];
	ierr = MPI_Send(buffer, 3*chunk, MPI_DOUBLE, 0, ibase, MPI_COMM_WORLD);
	if (ierr != 0) {
	  printf("Error in MPI_Send = %i\n",ierr);
	  ierr = MPI_Abort(MPI_COMM_WORLD, -1);
	  return -1;
	}

	/* update work counter */
	mycount += chunk;

      } /* if status.MPI_TAG */

    } /* end while loop */

    /* output work performed */
    printf(" proc %i computed %i iterations\n", myid, mycount);

  } /* end worker section */

  /* finalize MPI */
  ierr = MPI_Finalize();

} /* end main */
