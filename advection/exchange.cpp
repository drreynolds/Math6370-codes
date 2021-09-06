/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Implementation file for routines defined in exchange.h

// Inclusions
#include "advection_mpi.h"
#include "exchange.h"


// Initializes module data
void init_exchange(MPI_Comm com, int nxl, int nyl, 
		   int pcoords[2], exchange_data* exchange) {

  // declarations
  int nbcoords[2];
  int ierr;

  // set communicator
  exchange->comm = com;
    
  // set subdomain extents
  exchange->nxloc = nxl;
  exchange->nyloc = nyl;
    
  // allocate send buffers
  exchange->send_EW = new double[nyl];
  exchange->send_NS = new double[nxl];
    
  // determine process neighbors
  nbcoords[0] = pcoords[0]-1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(com, nbcoords, &(exchange->nbW));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"init_exchange error in MPI_Cart_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  nbcoords[0] = pcoords[0]+1;
  nbcoords[1] = pcoords[1];
  ierr = MPI_Cart_rank(com, nbcoords, &(exchange->nbE));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"init_exchange error in MPI_Cart_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]-1;
  ierr = MPI_Cart_rank(com, nbcoords, &(exchange->nbS));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"init_exchange error in MPI_Cart_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  nbcoords[0] = pcoords[0];
  nbcoords[1] = pcoords[1]+1;
  ierr = MPI_Cart_rank(com, nbcoords, &(exchange->nbN));
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"init_exchange error in MPI_Cart_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
}


// Frees exchange_data structure
void free_exchange(exchange_data* exchange) {
  delete[] exchange->send_EW;
  delete[] exchange->send_NS;
}


// Performs neighbor exchange to update v1 buffers
int do_exchange1(double* v1, double* v1W, double* v1S, 
		 exchange_data* exchange) {

  // check that module has been initialized
  if (exchange->nxloc < 0) {
    fprintf(stderr,"do_exchange1 error: module not initialized!\n");
    return 1;
  }

  // local variables
  MPI_Request reqW, reqS, sreqE, sreqN;
  MPI_Status status;
  int i, j, ierr;
  int nxl = exchange->nxloc;
  int nyl = exchange->nyloc;

  // open asynchronous receive channels to west, south neighbors
  ierr = MPI_Irecv(v1W, nyl, MPI_DOUBLE, exchange->nbW, 
		   1, exchange->comm, &reqW);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Irecv(v1S, nxl, MPI_DOUBLE, exchange->nbS, 
		   2, exchange->comm, &reqS);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }

  // pack send buffers
  for (j=0; j<nyl; j++)  
    exchange->send_EW[j] = v1[idx(nxl-1,j,nxl)];
  for (i=0; i<nxl; i++)  
    exchange->send_NS[i] = v1[idx(i,nyl-1,nxl)];

  // send data to neighbors
  ierr = MPI_Isend(exchange->send_EW, nyl, MPI_DOUBLE, 
		   exchange->nbE, 1, exchange->comm, &sreqE);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Isend = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Isend(exchange->send_NS, nxl, MPI_DOUBLE, 
		   exchange->nbN, 2, exchange->comm, &sreqN);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Isend = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  
  // wait for communications to finish
  ierr = MPI_Wait(&sreqE, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&sreqN, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&reqW, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&reqS, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange1 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }

  return 0;
}


// Performs neighbor exchange to update v2 and v3 buffers
int do_exchange2(double* v2, double* v3, double* v2E, 
		 double* v3N, exchange_data* exchange) {

  // check that module has been initialized
  if (exchange->nxloc < 0) {
    fprintf(stderr,"do_exchange2 error: module not initialized!\n");
    return 1;
  }
    
  // local variables
  MPI_Request reqE, reqN, sreqW, sreqS;
  MPI_Status status;
  int i, j, ierr;
  int nxl = exchange->nxloc;
  int nyl = exchange->nyloc;

  // open asynchronous receive channels to east, north
  ierr = MPI_Irecv(v2E, nyl, MPI_DOUBLE, exchange->nbE, 
		   1, exchange->comm, &reqE);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Irecv(v3N, nxl, MPI_DOUBLE, exchange->nbN, 
		   2, exchange->comm, &reqN);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Irecv = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }

  // pack send buffers
  for (j=0; j<nyl; j++)  
    exchange->send_EW[j] = v2[idx(0,j,nxl)];
  for (i=0; i<nxl; i++)  
    exchange->send_NS[i] = v3[idx(i,0,nxl)];

  // send data to neighbors
  ierr = MPI_Isend(exchange->send_EW, nyl, MPI_DOUBLE, 
		   exchange->nbW, 1, exchange->comm, &sreqW);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Isend = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Isend(exchange->send_NS, nxl, MPI_DOUBLE, 
		   exchange->nbS, 2, exchange->comm, &sreqS);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Isend = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
    
  // wait for communications to finish
  ierr = MPI_Wait(&sreqW, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&sreqS, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&reqE, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }
  ierr = MPI_Wait(&reqN, &status);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"do_exchange2 error in MPI_Wait = %i\n",ierr);
    MPI_Abort(exchange->comm, 1);
  }

  return 0;
}
