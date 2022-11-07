/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include "advection_mpi.hpp"


// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // intialize MPI and allocate parallel decomposition structure
  int ierr = MPI_Init(&argc, &argv);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Init");
  parallel_decomp p2d;

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(p2d.numprocs));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_size");

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(p2d.myid));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_rank");


  /* root process reads problem parameters from input file
     (should be in this order):
       nx - number of nodes in x-direction
       ny - number of nodes in y-direction
       nt - number of time steps
       tstop - final time (will stop at nt or stop, whichever is 1st)
       c - wave speed
       dtoutput - time frequency for outputting solutions
     and packs these into send buffers */
  double stime, ftime, iotime, tstop, c, dtoutput, dbuf[3];
  int nx, ny, nt, i, j, ibuf[3];
  if (p2d.myid == 0) {
    stime = MPI_Wtime();
    FILE *FID = fopen("input.txt","r");
    i = fscanf(FID," &inputs\n");
    i = fscanf(FID,"  nx = %i,\n", &nx);
    i = fscanf(FID,"  ny = %i,\n", &ny);
    i = fscanf(FID,"  nt = %i,\n", &nt);
    i = fscanf(FID,"  tstop = %lf,\n", &tstop);
    i = fscanf(FID,"  c = %lf,\n", &c);
    i = fscanf(FID,"  dtoutput = %lf,\n", &dtoutput);
    fclose(FID);
    printf("\nRunning wave problem:\n");
    printf("  nx = %i,  ny = %i\n", nx, ny);
    printf("  numprocs = %i\n", p2d.numprocs);
    printf("  nt = %i,  tstop = %g\n", nt, tstop);
    printf("  c = %g\n",c);
    printf("  dtoutput = %g\n",dtoutput);
    ftime = MPI_Wtime();
    iotime = ftime - stime;

    // fill Bcast buffers
    ibuf[0] = nx;
    ibuf[1] = ny;
    ibuf[2] = nt;
    dbuf[0] = tstop;
    dbuf[1] = c;
    dbuf[2] = dtoutput;
  }

  // broadcast data to all procs
  ierr = MPI_Bcast(&ibuf, 3, MPI_INT, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");
  ierr = MPI_Bcast(&dbuf, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");

  // unpack buffers
  nx = ibuf[0];
  ny = ibuf[1];
  nt = ibuf[2];
  tstop = dbuf[0];
  c = dbuf[1];
  dtoutput = dbuf[2];

  // set up parallel decomposition structure
  stime = MPI_Wtime();
  p2d.setup(nx, ny);

  // allocate arrays
  double *u  = new double[p2d.nxloc*p2d.nyloc];
  double *v1 = new double[p2d.nxloc*p2d.nyloc];
  double *v2 = new double[p2d.nxloc*p2d.nyloc];
  double *v3 = new double[p2d.nxloc*p2d.nyloc];

  // set grid spacing
  const double dx = 1.0/nx;
  const double dy = 1.0/ny;

  // set initial conditions
  initialize(u, v1, v2, v3, c, dx, dy, p2d);
  ftime = MPI_Wtime();
  double inittime = ftime - stime;

  // set initial time, output initial solution
  double t = 0.0;
  double toutput = 0.0;
  int noutput = 0;
  if (p2d.myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = MPI_Wtime();
  output(u, t, nx, ny, noutput, p2d);
  ftime = MPI_Wtime();
  iotime += (ftime - stime);

  // set time step
  const double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start time stepping
  double runtime = 0.0;
  double commtime = 0.0;
  for (int it=0; it<nt; it++) {

    // communicate v2 (EAST/WEST) & v3 (NORTH/SOUTH)
    stime = MPI_Wtime();
    ierr = p2d.Communication1(v2,v3);
    check_err(ierr, p2d.comm, "Communication1");
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    // first update v1 to get to half time step
    // get relevant values for this location
    stime = MPI_Wtime();
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	// access relevant components of v2 and v3
	double v2_E = (i == p2d.nxloc-1) ? p2d.v2recvE[j] : v2[idx(i+1,j,p2d.nxloc)];
        double v2_W = v2[idx(i,j,p2d.nxloc)];
        double v3_N = (j == p2d.nyloc-1) ? p2d.v3recvN[i] : v3[idx(i,j+1,p2d.nxloc)];
        double v3_S = v3[idx(i,j,p2d.nxloc)];

	// update v1
	v1[idx(i,j,p2d.nxloc)] += c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);

      } // for j
    } // for i
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    // communicate v1 (EAST/WEST) & v1 (NORTH/SOUTH)
    stime = MPI_Wtime();
    ierr = p2d.Communication2(v1);
    check_err(ierr, p2d.comm, "Communication2");
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    // next update v2 & v3 to get to full time step
    // get relevant values for this location
    stime = MPI_Wtime();
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	// access relevant components of v1
	double v1_W = (i == 0) ? p2d.v1recvW[j] : v1[idx(i-1,j,p2d.nxloc)];
        double v1_E = v1[idx(i,j,p2d.nxloc)];
        double v1_S = (j == 0) ? p2d.v1recvS[i] : v1[idx(i,j-1,p2d.nxloc)];
        double v1_N = v1[idx(i,j,p2d.nxloc)];

	// update v2 and v3
	v2[idx(i,j,p2d.nxloc)] += c*dt/dx*(v1_E - v1_W);
	v3[idx(i,j,p2d.nxloc)] += c*dt/dy*(v1_N - v1_S);

      } // for j
    } // for i

    // update solution for plotting
    for (j=0; j<p2d.nyloc; j++)
      for (i=0; i<p2d.nxloc; i++)
	u[idx(i,j,p2d.nxloc)] += dt*v1[idx(i,j,p2d.nxloc)];

    // update time
    t = t + dt;
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    // stop simulation if we've reached tstop
    if (t >= tstop)  break;

    // output solution periodically
    if ( t - (toutput + dtoutput) > -1.e-14 ) {
      stime = MPI_Wtime();
      toutput = t;
      noutput++;
      if (p2d.myid == 0)
        printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
      output(u, t, nx, ny, noutput, p2d);
      ftime = MPI_Wtime();
      iotime += ftime - stime;
    }

  } // for it


  // output final solution
  stime = MPI_Wtime();
  toutput = t;
  noutput++;
  if (p2d.myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, nt, t);
  output(u, t, nx, ny, noutput, p2d);
  ftime = MPI_Wtime();
  iotime += ftime - stime;

  // output times
  if (p2d.myid == 0) {
    printf(" total initialization time = %.16e\n", inittime);
    printf(" total input/output time   = %.16e\n", iotime);
    printf(" total communication time  = %.16e\n", commtime);
    printf(" total simulation time     = %.16e\n", runtime);
  }

  // free memory
  delete[] u;
  delete[] v1;
  delete[] v2;
  delete[] v3;
  p2d.free();

  // finalize MPI
  ierr = MPI_Finalize();

  return 0;
} // end main


//-------------------------------------------------------


// error checking routine for successful MPI calls
void check_err(const int ierr, MPI_Comm comm, const char* fname) {
  if (ierr != MPI_SUCCESS) {
    std::cerr << " error in " << fname << " = " << ierr << std::endl;
    MPI_Abort(comm, ierr);
  }
}


// setup routine for parallel decomposition structure
void parallel_decomp::setup(int nx, int ny) {

  // set up 2D Cartesian communicator
  this->periodic[0] = 1;
  this->periodic[1] = 1;
  this->pdims[0] = 0;
  this->pdims[1] = 0;
  int ierr = MPI_Dims_create(this->numprocs, 2, this->pdims);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Dims_create");

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, this->pdims, this->periodic, 0, &(this->comm));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Cart_create");

  ierr = MPI_Cart_coords(this->comm, this->myid, 2, this->pcoords);
  check_err(ierr, this->comm, "MPI_Cart_get");

  // determine local extents
  this->is = this->pcoords[0] * nx / this->pdims[0];
  this->ie = (this->pcoords[0] + 1) * nx / this->pdims[0] - 1;
  this->js = this->pcoords[1] * ny / this->pdims[1];
  this->je = (this->pcoords[1] + 1) * ny / this->pdims[1] - 1;

  // determine local array extents
  this->nxloc = this->ie - this->is + 1;
  this->nyloc = this->je - this->js + 1;

  // allocate send/recv arrays
  if (this->v2recvE == NULL)  this->v2recvE = new double[this->nyloc];
  if (this->v3recvN == NULL)  this->v3recvN = new double[this->nxloc];
  if (this->v2sendW == NULL)  this->v2sendW = new double[this->nyloc];
  if (this->v3sendS == NULL)  this->v3sendS = new double[this->nxloc];
  if (this->v1recvW == NULL)  this->v1recvW = new double[this->nyloc];
  if (this->v1recvS == NULL)  this->v1recvS = new double[this->nxloc];
  if (this->v1sendE == NULL)  this->v1sendE = new double[this->nyloc];
  if (this->v1sendN == NULL)  this->v1sendN = new double[this->nxloc];

  // determine process neighbors
  int nbcoords[2];
  nbcoords[0] = this->pcoords[0]-1;
  nbcoords[1] = this->pcoords[1];
  ierr = MPI_Cart_rank(this->comm, nbcoords, &(this->nbW));
  check_err(ierr, this->comm, "MPI_Cart_rank");

  nbcoords[0] = this->pcoords[0]+1;
  nbcoords[1] = this->pcoords[1];
  ierr = MPI_Cart_rank(this->comm, nbcoords, &(this->nbE));
  check_err(ierr, this->comm, "MPI_Cart_rank");

  nbcoords[0] = this->pcoords[0];
  nbcoords[1] = this->pcoords[1]-1;
  ierr = MPI_Cart_rank(this->comm, nbcoords, &(this->nbS));
  check_err(ierr, this->comm, "MPI_Cart_rank");

  nbcoords[0] = this->pcoords[0];
  nbcoords[1] = this->pcoords[1]+1;
  ierr = MPI_Cart_rank(this->comm, nbcoords, &(this->nbN));
  check_err(ierr, this->comm, "MPI_Cart_rank");
}


int parallel_decomp::Communication1(double *v2, double *v3) {

  // initialize send/recv buffers
  for (int j=0; j<nyloc; j++)  v2sendW[j] = v2[idx(0,j,nxloc)];
  for (int i=0; i<nxloc; i++)  v3sendS[i] = v3[idx(i,0,nxloc)];
  for (int j=0; j<nyloc; j++)  v2recvE[j] = 0.0;
  for (int i=0; i<nxloc; i++)  v3recvN[i] = 0.0;

  // open send/receive channels
  MPI_Request req[4];
  int ierr = MPI_Irecv(v2recvE, nyloc, MPI_DOUBLE, nbE, 2, comm, &(req[0]));
  check_err(ierr, comm, "MPI_Irecv");

  ierr = MPI_Irecv(v3recvN, nxloc, MPI_DOUBLE, nbN, 4, comm, &(req[1]));
  check_err(ierr, comm, "MPI_Irecv");

  ierr = MPI_Isend(v2sendW, nyloc, MPI_DOUBLE, nbW, 2, comm, &(req[2]));
  check_err(ierr, comm, "MPI_Isend");

  ierr = MPI_Isend(v3sendS, nxloc, MPI_DOUBLE, nbS, 4, comm, &(req[3]));
  check_err(ierr, comm, "MPI_Isend");

  // wait until all communications finish
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  check_err(ierr, comm, "MPI_Waitall");

  return 0;
}




int parallel_decomp::Communication2(double *v1) {

  // initialize send/recv buffers
  for (int j=0; j<nyloc; j++)  v1sendE[j] = v1[idx(nxloc-1,j,nxloc)];
  for (int i=0; i<nxloc; i++)  v1sendN[i] = v1[idx(i,nyloc-1,nxloc)];
  for (int j=0; j<nyloc; j++)  v1recvW[j] = 0.0;
  for (int i=0; i<nxloc; i++)  v1recvS[i] = 0.0;

  // open send/receive channels
  MPI_Request req[4];
  int ierr;
  ierr = MPI_Irecv(v1recvW, nyloc, MPI_DOUBLE, nbW, 2, comm, &(req[0]));
  check_err(ierr, comm, "MPI_Irecv");

  ierr = MPI_Irecv(v1recvS, nxloc, MPI_DOUBLE, nbS, 4, comm, &(req[1]));
  check_err(ierr, comm, "MPI_Irecv");

  ierr = MPI_Isend(v1sendE, nyloc, MPI_DOUBLE, nbE, 2, comm, &(req[2]));
  check_err(ierr, comm, "MPI_Isend");

  ierr = MPI_Isend(v1sendN, nxloc, MPI_DOUBLE, nbN, 4, comm, &(req[3]));
  check_err(ierr, comm, "MPI_Isend");

  // wait until all communications finish
  MPI_Status stat[4];
  ierr = MPI_Waitall(4, req, stat);
  check_err(ierr, comm, "MPI_Waitall");

  return 0;
}
