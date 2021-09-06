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
     Solves the linear system 
             a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k),
     using a parallelized direct tridiagonal solver. */
int tridiag_solve(double *a, double *b, double *c, double *u, 
		  double *r, int local_N, MPI_Comm comm) {

  // declarations
  int ierr, nprocs, my_id, k;
  double buff[3];
  double b_l, r_l, c_l, u_r;
  int msg_xch_fwd = 0;
  int msg_xch_bck = 1;
  MPI_Status status;

  // Get MPI parallelism information from comm
  ierr = MPI_Comm_size(comm, &nprocs);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"error in MPI_Comm_size = %i\n", ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Comm_rank(comm, &my_id);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr,"error in MPI_Comm_size = %i\n", ierr);
    MPI_Abort(comm, 1);
  }

  // FORWARD SWEEP
  //    Wait for left b,r,c values from left neighbor processor
  if (my_id > 0) {
    ierr = MPI_Recv(buff, 3, MPI_DOUBLE, my_id-1, msg_xch_fwd, comm, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Recv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
    b_l = buff[0];
    r_l = buff[1];
    c_l = buff[2];
  } else {
    b_l = 1.0;
    r_l = 0.0;
    c_l = 0.0;
  }
  
  //    Perform local sweep, updating diagonal and rhs
  b[0] -= c_l*a[0]/b_l;
  r[0] -= r_l*a[0]/b_l;
  for (k=1; k<local_N; k++) {
    b[k] -= c[k-1]*a[k]/b[k-1];
    r[k] -= r[k-1]*a[k]/b[k-1];
  }

  //    Send right-most b,r,c values to right neighbor processor
  if (my_id < nprocs-1) {
    buff[0] = b[local_N-1];
    buff[1] = r[local_N-1];
    buff[2] = c[local_N-1];
    ierr = MPI_Send(buff, 3, MPI_DOUBLE, my_id+1, msg_xch_fwd, comm);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Send = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }


  // BACKWARD SWEEP
  //    Wait for right boundary value from right neighbor processor
  if (my_id < nprocs-1) {
    ierr = MPI_Recv(&u_r, 1, MPI_DOUBLE, my_id+1, msg_xch_bck, comm, &status);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Recv = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  } else {
    u_r = 0.0;
  }

  //    Perform local sweep, updating solution
  u[local_N-1] = (r[local_N-1] - c[local_N-1]*u_r)/b[local_N-1];
  for (k=local_N-2; k>=0; k--)
    u[k] = (r[k]-c[k]*u[k+1])/b[k];

  //    Send left-most value to left neighbor processor
  if (my_id > 0) {
    u_r = u[0];
    ierr = MPI_Send(&u_r, 1, MPI_DOUBLE, my_id-1, msg_xch_bck, comm);
    if (ierr != 0) {
      fprintf(stderr,"error in MPI_Send = %i\n",ierr);
      MPI_Abort(comm, 1);
    }
  }


  // Wait until all processes have caught up
  ierr = MPI_Barrier(comm);
  if (ierr != 0) {
    fprintf(stderr,"error in MPI_Barrier = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  
  // finished
  return 0;

} // end tridiag_solve
