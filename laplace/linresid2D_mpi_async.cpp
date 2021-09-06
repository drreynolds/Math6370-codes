/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   28 February 2013 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "laplace2d.h"


/* Description:
      calculates the 2D linear residual 
               res = L*u - f
      and its L2-norm. */
int linresid2D(double *u, double *f, double *res, double *norm2, 
	       int locN, int locM, double dx, double dy, MPI_Comm comm) {

  // declarations
  int ierr=0, my_id, i, j, nbW, nbE, nbS, nbN;
  int pdims[2], pcoords[2], periods[2], nbcoords[2];
  double norm;
  double *Esend, *Wsend, *Nsend, *Ssend, *Erecv, *Wrecv, *Nrecv, *Srecv;
  MPI_Request reqE, reqN, reqW, reqS, sreqE, sreqN, sreqW, sreqS;
  MPI_Status status;

  // get MPI parallelism information from comm
  ierr = MPI_Cart_get(comm, 2, pdims, periods, pcoords);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_get = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Comm_rank(comm, &my_id);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // initialize send/recv buffers
  Esend = new double[locM];
  Wsend = new double[locM];
  Nsend = new double[locN];
  Ssend = new double[locN];
  Erecv = new double[locM];
  Wrecv = new double[locM];
  Nrecv = new double[locN];
  Srecv = new double[locN];
  for (j=0; j<locM; j++)  Esend[j] = u[idx(locN-1,j,locN)];
  for (j=0; j<locM; j++)  Wsend[j] = u[idx(0,j,locN)];
  for (i=0; i<locN; i++)  Nsend[i] = u[idx(i,locM-1,locN)];
  for (i=0; i<locN; i++)  Ssend[i] = u[idx(i,0,locN)];
  for (j=0; j<locM; j++)  Erecv[j] = 0.0;
  for (j=0; j<locM; j++)  Wrecv[j] = 0.0;
  for (i=0; i<locN; i++)  Nrecv[i] = 0.0;
  for (i=0; i<locN; i++)  Srecv[i] = 0.0;

  // determine process 'neighbors'
  if (pcoords[0] > 0) {
    nbcoords[0] = pcoords[0]-1;
    nbcoords[1] = pcoords[1];
    ierr = MPI_Cart_rank(comm, nbcoords, &nbW);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[0] < pdims[0]-1) {
    nbcoords[0] = pcoords[0]+1;
    nbcoords[1] = pcoords[1];
    ierr = MPI_Cart_rank(comm, nbcoords, &nbE);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] > 0) {
    nbcoords[0] = pcoords[0];
    nbcoords[1] = pcoords[1]-1;
    ierr = MPI_Cart_rank(comm, nbcoords, &nbS);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] < pdims[1]-1) {
    nbcoords[0] = pcoords[0];
    nbcoords[1] = pcoords[1]+1;
    ierr = MPI_Cart_rank(comm, nbcoords, &nbN);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Cart_rank = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // phase 1: open receive channels for neighbor values
  if (pcoords[0] > 0) {
    ierr = MPI_Irecv(Wrecv, locM, MPI_DOUBLE, nbW, 1, comm, &reqW);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[0] < pdims[0]-1) {
    ierr = MPI_Irecv(Erecv, locM, MPI_DOUBLE, nbE, 2, comm, &reqE);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] > 0) {
    ierr = MPI_Irecv(Srecv, locN, MPI_DOUBLE, nbS, 3, comm, &reqS);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] < pdims[1]-1) {
    ierr = MPI_Irecv(Nrecv, locN, MPI_DOUBLE, nbN, 4, comm, &reqN);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // phase 2: send boundary values to neighbors
  if (pcoords[0] < pdims[0]-1) {
    ierr = MPI_Isend(Esend, locM, MPI_DOUBLE, nbE, 1, comm, &sreqE);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[0] > 0) {
    ierr = MPI_Isend(Wsend, locM, MPI_DOUBLE, nbW, 2, comm, &sreqW);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] < pdims[1]-1) {
    ierr = MPI_Isend(Nsend, locN, MPI_DOUBLE, nbN, 3, comm, &sreqN);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] > 0) {
    ierr = MPI_Isend(Ssend, locN, MPI_DOUBLE, nbS, 4, comm, &sreqS);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }


  // phase 3: compute linear residual in interior of subdomain
  norm = 0.0;
  for (i=1; i<locN-1; i++) {
    for (j=1; j<locM-1; j++) {
      res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
	+ (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
	+ (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
      norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
    }
  }


  // phase 4: wait until receives finish
  if (pcoords[0] > 0) {
    ierr = MPI_Wait(&reqW, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[0] < pdims[0]-1) {
    ierr = MPI_Wait(&reqE, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] > 0) {
    ierr = MPI_Wait(&reqS, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] < pdims[1]-1) {
    ierr = MPI_Wait(&reqN, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }


  // phase 5: compute residual along edges

  //   compute linear residual at South edge of local subdomain
  j=0;
  for (i=1; i<locN-1; i++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
    norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  //   compute linear residual at West edge of local subdomain
  i=0;
  for (j=1; j<locM-1; j++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)] 
      + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
    norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }
    
  //   compute linear residual at East edge of local subdomain
  i=locN-1;
  for (j=1; j<locM-1; j++) {
    res[idx(i,j,locN)] = -f[idx(i,j,locN)] 
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy ;
    norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  //   compute linear residual at North edge of local subdomain
  j=locM-1;
  for (i=1; i<locN-1; i++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)] 
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
    norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  // phase 6: compute linear residual at corners of local subdomain
  i=0; j=0;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)] 
    + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
    + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
  norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=locN-1; j=0;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
    + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
  norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=0; j=locM-1;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
    + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
  norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=locN-1; j=locM-1;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
    + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
  norm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];


  // phase 7: combine local sums into global L2-norm
  ierr = MPI_Allreduce(&norm, norm2, 1, MPI_DOUBLE, MPI_SUM, comm);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Allreduce = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  *norm2 = sqrt(*norm2);

  // phase 8: wait until sends finish
  if (pcoords[0] > 0) {
    ierr = MPI_Wait(&sreqW, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[0] < pdims[0]-1) {
    ierr = MPI_Wait(&sreqE, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] > 0) {
    ierr = MPI_Wait(&sreqS, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }
  if (pcoords[1] < pdims[1]-1) {
    ierr = MPI_Wait(&sreqN, &status);
    if (ierr != MPI_SUCCESS) {
      fprintf(stderr," error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }


  // delete temporary buffers
  delete[] Esend;
  delete[] Wsend;
  delete[] Nsend;
  delete[] Ssend;
  delete[] Erecv;
  delete[] Wrecv;
  delete[] Nrecv;
  delete[] Srecv;

  return ierr;

} // end linresid

