/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

/* Example routine using the basic 6 MPI functions */
int main(int argc, char* argv[]) {

  /* local variables */
  int numprocs, myid, number, p, tag, sender;
  MPI_Status status;

  /* intialize MPI */
  if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Init\n");
    return 1;
  }
  if (MPI_Comm_size(MPI_COMM_WORLD, &numprocs) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Comm_size\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
  if (MPI_Comm_rank(MPI_COMM_WORLD, &myid) != MPI_SUCCESS) {
    fprintf(stderr,"Error in MPI_Comm_rank\n");
    return 1;
  }

  /* everyone send their (ID+1) to the root processor (except for root)
     - set the tag to mirror the sending processor id */
  if (myid != 0) {
    number = myid + 1;
    if (MPI_Send(&number, 1, MPI_INT, 0, myid,
		 MPI_COMM_WORLD) != MPI_SUCCESS) {
      fprintf(stderr,"Error in MPI_Send\n");
      return 1;
    }
  }

  /* the root node receives these (in order) and outputs each to screen */
  if (myid == 0) {

    /* loop over all other processors */
    for (p=1; p<numprocs; p++) {

      /* receive the number from this processor */
      if (MPI_Recv(&number, 1, MPI_INT, p, MPI_ANY_TAG,
		   MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
	fprintf(stderr,"Error in MPI_Recv\n");
	return 1;
      }

      /* get some information about this message */
      tag    = status.MPI_TAG;
      sender = status.MPI_SOURCE;

      /* output the information and the data to screen */
      printf("received value %i from processor %i, tag = %i, sender = %i\n",
	     number, p, tag, sender);

    } /* for p */
  } /* if myid */

  /* finalize MPI */
  MPI_Finalize();

  return 0;
} /* end main */
