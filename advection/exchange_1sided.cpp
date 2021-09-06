/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */


// Inclusions 
#include <stdlib.h>     // new, delete
#include <math.h>       // exp(), pow()
#include "advection_1sided.h"


// creates RMA windows into recv buffers for subsequent calls to MPI_Rput operations
int setup_windows(int nxl, int nyl, double* v1_recvW, 
		  double* v1_recvS, double* v2_recvE, 
		  double* v3_recvN, MPI_Comm comm, 
		  MPI_Win& win_v1W, MPI_Win& win_v1S, 
		  MPI_Win& win_v2E, MPI_Win& win_v3N,
		  MPI_Group& groupW, MPI_Group& groupE, 
		  MPI_Group& groupS, MPI_Group& groupN) {

  // local data
  int ierr, nbN, nbS, nbE, nbW;
  int pdims[2], periods[2], pcoords[2], nbcoords[2];
  MPI_Group comm_group;

  // figure out where we are, who our neighbors are
  ierr = MPI_Cart_get(comm, 2, pdims, periods, pcoords);
  check_err(ierr, comm, "MPI_Cart_get");

  nbcoords[0] = pcoords[0]-1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(comm, nbcoords, &nbW);
  check_err(ierr, comm, "MPI_Cart_rank");

  nbcoords[0] = pcoords[0]+1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(comm, nbcoords, &nbE);
  check_err(ierr, comm, "MPI_Cart_rank");

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]-1;
  ierr = MPI_Cart_rank(comm, nbcoords, &nbS);
  check_err(ierr, comm, "MPI_Cart_rank");

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]+1;
  ierr = MPI_Cart_rank(comm, nbcoords, &nbN);
  check_err(ierr, comm, "MPI_Cart_rank");

  // create MPI group from full communicator
  ierr = MPI_Comm_group(comm, &comm_group);
  check_err(ierr, comm, "MPI_Comm_group");

  // create nearest-neighbor MPI send/receive groups
  ierr = MPI_Group_incl(comm_group, 1, &nbW, &groupW);
  check_err(ierr, comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &nbE, &groupE);
  check_err(ierr, comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &nbS, &groupS);
  check_err(ierr, comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &nbN, &groupN);
  check_err(ierr, comm, "MPI_Group_incl");

  // create RMA window into v1_recvW
  ierr = MPI_Win_create(v1_recvW, nyl*sizeof(double), sizeof(double), 
			MPI_INFO_NULL, comm, &win_v1W);
  check_err(ierr, comm, "MPI_Win_create");

  // create RMA window into v1_recvS
  ierr = MPI_Win_create(v1_recvS, nxl*sizeof(double), sizeof(double), 
			MPI_INFO_NULL, comm, &win_v1S);
  check_err(ierr, comm, "MPI_Win_create");

  // create RMA window into v2_recvE
  ierr = MPI_Win_create(v2_recvE, nyl*sizeof(double), sizeof(double), 
			MPI_INFO_NULL, comm, &win_v2E);
  check_err(ierr, comm, "MPI_Win_create");

  // create RMA window into v3_recvN
  ierr = MPI_Win_create(v3_recvN, nxl*sizeof(double), sizeof(double), 
			MPI_INFO_NULL, comm, &win_v3N);
  check_err(ierr, comm, "MPI_Win_create");

  // free general communicator group
  ierr = MPI_Group_free(&comm_group);
  check_err(ierr, comm, "MPI_Group_free");

  // return success
  return 0;
}


// fill neighbors' v2_recvE and v3_recvN buffers using MPI_Put calls
int put_v2v3(double* v2, double* v3, double* v2_sendW, 
	     double* v3_sendS, int nxl, int nyl, MPI_Comm comm, 
	     MPI_Win& win_v2E, MPI_Win& win_v3N) {

  // local data
  int i, j, ierr, nbW, nbS;
  int pdims[2], periods[2], pcoords[2], nbcoords[2];

  // fill send buffers
  for (j=0; j<nyl; j++) 
    v2_sendW[j] = v2[idx(0,j,nxl)];
  for (i=0; i<nxl; i++) 
    v3_sendS[i] = v3[idx(i,0,nxl)];
  
  // figure out where we are, who W and S neighbors are
  ierr = MPI_Cart_get(comm, 2, pdims, periods, pcoords);
  check_err(ierr, comm, "MPI_Cart_get");

  nbcoords[0] = pcoords[0]-1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(comm, nbcoords, &nbW);
  check_err(ierr, comm, "MPI_Cart_rank");

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]-1;
  ierr = MPI_Cart_rank(comm, nbcoords, &nbS);
  check_err(ierr, comm, "MPI_Cart_rank");
  
  // start MPI_Put calls to put v2_sendW and v3_sendS into neighbor buffers
  ierr = MPI_Put(v2_sendW, nyl, MPI_DOUBLE, nbW, 0, nyl, MPI_DOUBLE, win_v2E);
  check_err(ierr, comm, "MPI_Put");

  ierr = MPI_Put(v3_sendS, nxl, MPI_DOUBLE, nbS, 0, nxl, MPI_DOUBLE, win_v3N);
  check_err(ierr, comm, "MPI_Put");

  return 0;
}

// fill neighbors' v1_recvW and v1_recvS buffers using MPI_Put calls
int put_v1(double* v1, double* v1_sendE, double* v1_sendN, int nxl, 
	   int nyl, MPI_Comm comm, MPI_Win& win_v1W, MPI_Win& win_v1S) {

  // local data
  int i, j, ierr, nbE, nbN;
  int pdims[2], periods[2], pcoords[2], nbcoords[2];

  // fill send buffers
  for (j=0; j<nyl; j++) 
    v1_sendE[j] = v1[idx(nxl-1,j,nxl)];
  for (i=0; i<nxl; i++) 
    v1_sendN[i] = v1[idx(i,nyl-1,nxl)];
  
  // figure out where we are, E and N neighbors are
  ierr = MPI_Cart_get(comm, 2, pdims, periods, pcoords);
  check_err(ierr, comm, "MPI_Cart_get");

  nbcoords[0] = pcoords[0]+1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(comm, nbcoords, &nbE);
  check_err(ierr, comm, "MPI_Cart_rank");

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]+1;
  ierr = MPI_Cart_rank(comm, nbcoords, &nbN);
  check_err(ierr, comm, "MPI_Cart_rank");

  // start MPI_Put calls to send v1_sendE and v1_sendN
  ierr = MPI_Put(v1_sendE, nyl, MPI_DOUBLE, nbE, 0, nyl, MPI_DOUBLE, win_v1W);
  check_err(ierr, comm, "MPI_Put");

  ierr = MPI_Put(v1_sendN, nxl, MPI_DOUBLE, nbN, 0, nxl, MPI_DOUBLE, win_v1S);
  check_err(ierr, comm, "MPI_Put");

  return 0;
}

