/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   28 February 2013 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"


/* Description: 
      calculates the linear residual and returns its averaged 2-norm (WRMS) */
double linresid(double *a, double *b, double *c, double *u, double *r, 
		double *res, int local_N, int global_N, MPI_Comm comm) {

  // declarations
  int ierr, nprocs, my_id, k;
  double u_l, u_r, s_l, s_r, tmp, norm2;
  MPI_Request id_recv_l, id_recv_r, id_send_l, id_send_r;
  MPI_Status status;

  // Get MPI parallelism information from comm
  ierr = MPI_Comm_size(comm, &nprocs);
  if (ierr != 0) {
    fprintf(stderr,"error in MPI_Comm_size = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Comm_rank(comm, &my_id);
  if (ierr != 0) {
    fprintf(stderr,"error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // open asynchronous receive channels for neighbor boundary values
  u_l = 0.0;
  if (my_id > 0) {
    ierr = MPI_Irecv(&u_l, 1, MPI_DOUBLE, my_id-1, 100, comm, &id_recv_l);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  } 
  u_r = 0.0;
  if (my_id < nprocs-1) {
    ierr = MPI_Irecv(&u_r, 1, MPI_DOUBLE, my_id+1, 101, comm, &id_recv_r);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Irecv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // send boundary values to neighbor processes
  s_r = u[local_N-1];
  if (my_id < nprocs-1) {
    ierr = MPI_Isend(&s_r, 1, MPI_DOUBLE, my_id+1, 100, comm, &id_send_r);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  } 
  s_l = u[0];
  if (my_id > 0) {
    ierr = MPI_Isend(&s_l, 1, MPI_DOUBLE, my_id-1, 101, comm, &id_send_l);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Isend = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // compute linear residual in interior of subdomain
  norm2 = 0.0;
  for (k=1; k<local_N-1; k++) {
    res[k] = a[k]*u[k-1] + b[k]*u[k] + c[k]*u[k+1] - r[k];
    norm2 += res[k]*res[k];
  }

  // wait for left boundary value to arrive before using
  if (my_id > 0) {
    ierr = MPI_Wait(&id_recv_l, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // wait for right boundary value to arrive before using
  if (my_id < nprocs-1) {
    ierr = MPI_Wait(&id_recv_r, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // compute linear residual at left of subdomain
  res[0] = a[0]*u_l + b[0]*u[0] + c[0]*u[1] - r[0];
  norm2 += res[0]*res[0];

  // compute linear residual at right end of subdomain
  k = local_N-1;
  res[k] = a[k]*u[k-1] + b[k]*u[k]+ c[k]*u_r - r[k];
  norm2 += res[k]*res[k];

  // combine local 2-norms into global averaged 2-norm
  tmp = 0.0;
  ierr = MPI_Allreduce(&norm2, &tmp, 1, MPI_DOUBLE, MPI_SUM, comm);
  if (ierr != 0) {
    fprintf(stderr,"error in MPI_Allreduce = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  norm2 = sqrt(tmp/global_N);


  // wait for sends to finish
  if (my_id > 0) {
    ierr = MPI_Wait(&id_send_l, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  // wait for right boundary value to arrive before using
  if (my_id < nprocs-1) {
    ierr = MPI_Wait(&id_send_r, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Wait = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }

  return norm2;

} // end linresid
