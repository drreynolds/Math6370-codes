/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "laplace2d.hpp"


/*   We set up the linear residual
             L*u - f,
     where L is a standard 2D Laplace operator, and f and u are
     given vectors.  These are decomposed using a 2D parallel
     domain decomposition strategy.  We then call the routine
     linresid2D to compute the linear residual above. */
int main(int argc, char* argv[]) {

  // intialize MPI
  int ierr, numprocs, myid;
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_size = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  /* read problem parameters from input file (should be in this order):
         N - number of nodes in x-direction
	 M - number of nodes in y-direction
	px - number of procs in x-direction
	py - number of procs in y-direction */
  int N, M, px, py, buf[4];
  FILE *FID;
  if (myid == 0) {
    FID=fopen("input.txt","r");
    ierr = fscanf(FID," &inputs\n");
    ierr = fscanf(FID,"  N = %i,\n", &N);
    ierr = fscanf(FID,"  M = %i,\n", &M);
    ierr = fscanf(FID,"  px = %i,\n", &px);
    ierr = fscanf(FID,"  py = %i,\n", &py);
    fclose(FID);

    // fill send buffer
    buf[0] = N;
    buf[1] = M;
    buf[2] = px;
    buf[3] = py;
  }

  // broadcast data to all procs
  ierr = MPI_Bcast(&buf, 4, MPI_INT, 0, MPI_COMM_WORLD);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Bcast = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // extract values from buffer
  N  = buf[0];
  M  = buf[1];
  px = buf[2];
  py = buf[3];

  // check that the processor layout and communicator size agree
  if (numprocs != px*py) {
    if (myid == 0) {
      fprintf(stderr," error: incorrect processor layout\n");
      fprintf(stderr,"     numprocs = %i\n",numprocs);
      fprintf(stderr,"     px = %i,  py = %i\n",px,py);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // set up Cartesian communicator, comm
  int pdims[2], periods[2];
  MPI_Comm comm;
  pdims[0] = px;
  pdims[1] = py;
  periods[0] = 0;
  periods[1] = 0;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, 0, &comm);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_create = %i\n",ierr);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // get this processor's new id and location from comm
  int pcoords[2];
  ierr = MPI_Comm_rank(comm, &myid);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Comm_rank = %i\n",ierr);
    MPI_Abort(comm, 1);
  }
  ierr = MPI_Cart_coords(comm, myid, 2, pcoords);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Cart_coords = %i\n",ierr);
    MPI_Abort(comm, 1);
  }

  // set local mesh sizes (last process in each dir takes up remainder)
  const double dx = 1.0/(N-1);
  const double dy = 1.0/(M-1);
  int locN = N/px;
  int locM = M/py;
  const double xl = dx*locN*pcoords[0];
  const double yl = dy*locM*pcoords[1];
  if (pcoords[0] == px-1)
    if (locN*px != N)
      locN = N - locN*(px-1);
  if (pcoords[1] == py-1)
    if (locM*py != M)
      locM = M - locM*(py-1);

  // root node outputs some information to screen
  if (myid == 0) {
    printf("Laplace2D driver with %i processes\n", numprocs);
    printf("    global problem size = %i by %i\n", N, M);
    printf("    local problem sizes = %i by %i\n", locN, locM);
  }

  // allocate vector memory
  double *u   = new double[locN*locM];
  double *f   = new double[locN*locM];
  double *res = new double[locN*locM];

  // start timer
  double stime = MPI_Wtime();

  // set up u and f
  int i, j;
  for (j=0; j<locM; j++) {
    double y = yl + j*dy;
    for (i=0; i<locN; i++) {
      double x = xl + i*dx;
      u[idx(i,j,locN)] = exp(-20.0*((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)));
      f[idx(i,j,locN)] = 1.0;
    }
  }

  // adjust u at domain boundary
  if (pcoords[0] == 0)
    for (j=0; j<locM; j++)
      u[idx(0,j,locN)] = 0.0;
  if (pcoords[0] == px-1)
    for (j=0; j<locM; j++)
      u[idx(locN-1,j,locN)] = 0.0;
  if (pcoords[1] == 0)
    for (i=0; i<locN; i++)
      u[idx(i,0,locN)] = 0.0;
  if (pcoords[1] == py-1)
    for (i=0; i<locN; i++)
      u[idx(i,locM-1,locN)] = 0.0;

  // wait until all procs have caught up
  ierr = MPI_Barrier(comm);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Barrier = %i\n",ierr);
    ierr = MPI_Abort(comm, 1);
  }

  /* compute linear residual */
  double res2norm;
  ierr = linresid2D(u, f, res, res2norm, locN, locM, dx, dy, comm);
  if (ierr != MPI_SUCCESS) {
    printf(" error in linresid2D = %i\n",ierr);
    ierr = MPI_Abort(comm, 1);
  }
  if (myid == 0)
    printf(" residual: ||L*u-f||_2 = %g\n",res2norm);

  // stop timer
  double ftime = MPI_Wtime();
  if (myid == 0)
    printf(" runtime = %g\n",ftime-stime);

  // output residual to file(s)
  if (myid == 0) {
    FID = fopen("resid.txt","w");
    fprintf(FID, "%i %i %i %i\n", N, M, px, py);
    fclose(FID);
  }
  char outname[100];
  sprintf(outname, "resid.%03i", myid);
  FID = fopen(outname,"w");
  fprintf(FID, "%i\n", locN);
  fprintf(FID, "%i\n", locM);
  fprintf(FID, "%i\n", pcoords[0]);
  fprintf(FID, "%i\n", pcoords[1]);
  for (j=0; j<locM; j++)
    for (i=0; i<locN; i++)
      fprintf(FID, "%.16e\n",res[idx(i,j,locN)]);
  fclose(FID);

  // free memory
  delete[] u;
  delete[] f;
  delete[] res;

  // finalize MPI
  MPI_Finalize();

} // end main