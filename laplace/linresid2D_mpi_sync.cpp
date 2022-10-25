/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

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
int linresid2D(double *u, double *f, double *res, double &norm,
	       int locN, int locM, double dx, double dy, MPI_Comm comm) {

  // get MPI parallelism information from comm
  int ierr, pdims[2], pcoords[2], periods[2];
  ierr = MPI_Cart_get(comm, 2, pdims, periods, pcoords);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_get = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  /* allocate send/recv buffers */
  double *Esend = new double[locM];
  double *Wsend = new double[locM];
  double *Nsend = new double[locN];
  double *Ssend = new double[locN];
  double *Erecv = new double[locM];
  double *Wrecv = new double[locM];
  double *Nrecv = new double[locN];
  double *Srecv = new double[locN];

  /* initialize send/recv buffers */
  int i, j;
  for (j=0; j<locM; j++)  Esend[j] = u[idx(locN-1,j,locN)];
  for (j=0; j<locM; j++)  Wsend[j] = u[idx(0,j,locN)];
  for (i=0; i<locN; i++)  Nsend[i] = u[idx(i,locM-1,locN)];
  for (i=0; i<locN; i++)  Ssend[i] = u[idx(i,0,locN)];
  for (j=0; j<locM; j++)  Erecv[j] = 0.0;
  for (j=0; j<locM; j++)  Wrecv[j] = 0.0;
  for (i=0; i<locN; i++)  Nrecv[i] = 0.0;
  for (i=0; i<locN; i++)  Srecv[i] = 0.0;

  // determine process 'neighbors'
  int nbW, nbE, nbS, nbN, nbcoords[2];
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

  // phase 1: even pcoords[0] send East, odd recv West
  MPI_Status status;
  if (pcoords[0]%2 == 0) {
    if (pcoords[0] < pdims[0]-1) {
      ierr = MPI_Send(Esend, locM, MPI_DOUBLE, nbE, 100, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[0] > 0) {
      ierr = MPI_Recv(Wrecv, locM, MPI_DOUBLE, nbW, 100, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 2: even pcoords[0] recv East, odd send West
  if (pcoords[0]%2 == 0) {
    if (pcoords[0] < pdims[0]-1) {
      ierr = MPI_Recv(Erecv, locM, MPI_DOUBLE, nbE, 101, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[0] > 0) {
      ierr = MPI_Send(Wsend, locM, MPI_DOUBLE, nbW, 101, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 3: odd pcoords[0] send East, even recv West
  if (pcoords[0]%2 == 1) {
    if (pcoords[0] < pdims[0]-1) {
      ierr = MPI_Send(Esend, locM, MPI_DOUBLE, nbE, 102, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[0] > 0) {
      ierr = MPI_Recv(Wrecv, locM, MPI_DOUBLE, nbW, 102, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 4: odd pcoords[0] recv East, even send West
  if (pcoords[0]%2 == 1) {
    if (pcoords[0] < pdims[0]-1) {
      ierr = MPI_Recv(Erecv, locM, MPI_DOUBLE, nbE, 103, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[0] > 0) {
      ierr = MPI_Send(Wsend, locM, MPI_DOUBLE, nbW, 103, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 5: even pcoords[1] send North, odd recv South
  if (pcoords[1]%2 == 0) {
    if (pcoords[1] < pdims[1]-1) {
      ierr = MPI_Send(Nsend, locN, MPI_DOUBLE, nbN, 104, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[1] > 0) {
      ierr = MPI_Recv(Srecv, locN, MPI_DOUBLE, nbS, 104, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 6: even pcoords[1] recv North, odd send South
  if (pcoords[1]%2 == 0) {
    if (pcoords[1] < pdims[1]-1) {
      ierr = MPI_Recv(Nrecv, locN, MPI_DOUBLE, nbN, 105, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[1] > 0) {
      ierr = MPI_Send(Ssend, locN, MPI_DOUBLE, nbS, 105, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 7: odd pcoords[1] send North, even recv South
  if (pcoords[1]%2 == 1) {
    if (pcoords[1] < pdims[1]-1) {
      ierr = MPI_Send(Nsend, locN, MPI_DOUBLE, nbN, 106, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[1] > 0) {
      ierr = MPI_Recv(Srecv, locN, MPI_DOUBLE, nbS, 106, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // phase 8: odd pcoords[1] recv North, even send South
  if (pcoords[1]%2 == 1) {
    if (pcoords[1] < pdims[1]-1) {
      ierr = MPI_Recv(Nrecv, locN, MPI_DOUBLE, nbN, 107, comm, &status);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Recv = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  } else {
    if (pcoords[1] > 0) {
      ierr = MPI_Send(Ssend, locN, MPI_DOUBLE, nbS, 107, comm);
      if (ierr != MPI_SUCCESS) {
	fprintf(stderr," error in MPI_Send = %i\n",ierr);
	MPI_Abort(comm, 1);
      }
    }
  }

  // compute linear residual in interior of subdomain
  double mynorm = 0.0;
  for (i=1; i<locN-1; i++) {
    for (j=1; j<locM-1; j++) {
      res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
	+ (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
	+ (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
      mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
    }
  }

  // compute linear residual at South edge of local subdomain
  j=0;
  for (i=1; i<locN-1; i++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
    mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  // compute linear residual at West edge of local subdomain
  i=0;
  for (j=1; j<locM-1; j++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
      + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
    mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  // compute linear residual at East edge of local subdomain
  i=locN-1;
  for (j=1; j<locM-1; j++) {
    res[idx(i,j,locN)] = -f[idx(i,j,locN)]
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy ;
    mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  // compute linear residual at North edge of local subdomain
  j=locM-1;
  for (i=1; i<locN-1; i++) {
    res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
      + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
      + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
    mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];
  }

  // compute linear residual at corners of local subdomain
  i=0; j=0;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
    + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
  mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=locN-1; j=0;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
    + (Srecv[i] - 2.0*u[idx(i,j,locN)] + u[idx(i,j+1,locN)])/dy/dy;
  mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=0; j=locM-1;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (Wrecv[j] - 2.0*u[idx(i,j,locN)] + u[idx(i+1,j,locN)])/dx/dx
    + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
  mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];

  i=locN-1; j=locM-1;
  res[idx(i,j,locN)] =  -f[idx(i,j,locN)]
    + (u[idx(i-1,j,locN)] - 2.0*u[idx(i,j,locN)] + Erecv[j])/dx/dx
    + (u[idx(i,j-1,locN)] - 2.0*u[idx(i,j,locN)] + Nrecv[i])/dy/dy;
  mynorm += dx*dy*res[idx(i,j,locN)]*res[idx(i,j,locN)];


  // combine local sums into global L2-norm
  ierr = MPI_Allreduce(&mynorm, &norm, 1, MPI_DOUBLE, MPI_SUM, comm);
  if (ierr != MPI_SUCCESS) {
    printf(" error in MPI_Allreduce = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  norm = sqrt(norm);


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
