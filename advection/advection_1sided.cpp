/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cmath>
#include "advection_1sided.hpp"


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
    FILE* FID = fopen("input.txt","r");
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
    printf("  nt = %i,  tstop = %g\n", nt, tstop);
    printf("  numprocs = %i\n", p2d.numprocs);
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
  iotime += ftime - stime;

  // set time step
  const double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start time stepping
  double runtime = 0.0;
  double commtime = 0.0;
  for (int it=0; it<nt; it++) {

    // perform v2v3 communication
    stime = MPI_Wtime();    // start timer
    p2d.start_v2v3();       // expose receive buffers, start RMA access epochs
    p2d.put_v2v3(v2, v3);   // put v2 and v3 values into neighbor buffers
    p2d.finish_v2v3();      // stop RMA access epochs, wait for our buffers to fill
    ftime = MPI_Wtime();    // stop/increment timer
    commtime += ftime - stime;

    // first update v1 to get to half time step
    // get relevant values for this location
    stime = MPI_Wtime();
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	// access relevant components of v2 and v3
	double v2_E = (i == p2d.nxloc-1) ? p2d.v2_recvE[j] : v2[idx(i+1,j,p2d.nxloc)];
        double v2_W = v2[idx(i,j,p2d.nxloc)];
        double v3_N = (j == p2d.nyloc-1) ? p2d.v3_recvN[i] : v3[idx(i,j+1,p2d.nxloc)];
        double v3_S = v3[idx(i,j,p2d.nxloc)];

	// update v1
	v1[idx(i,j,p2d.nxloc)] += c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);

      } // for j
    } // for i
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    // perform v1 communication
    stime = MPI_Wtime();    // start timer
    p2d.start_v1();         // expose receive buffers, start RMA access epochs
    p2d.put_v1(v1);         // put v1 values into neighbor buffers
    p2d.finish_v1();        // stop RMA access epochs, wait for our buffers to fill
    ftime = MPI_Wtime();    // stop/increment timer
    commtime += ftime - stime;

    // next update v2 & v3 to get to full time step
    // get relevant values for this location
    stime = MPI_Wtime();
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	// access relevant components of v1
	double v1_W = (i == 0) ? p2d.v1_recvW[j] : v1[idx(i-1,j,p2d.nxloc)];
        double v1_E = v1[idx(i,j,p2d.nxloc)];
        double v1_S = (j == 0) ? p2d.v1_recvS[i] : v1[idx(i,j-1,p2d.nxloc)];
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
  MPI_Finalize();

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
  if (this->v2_recvE == NULL)  this->v2_recvE = new double[this->nyloc];
  if (this->v3_recvN == NULL)  this->v3_recvN = new double[this->nxloc];
  if (this->v2_sendW == NULL)  this->v2_sendW = new double[this->nyloc];
  if (this->v3_sendS == NULL)  this->v3_sendS = new double[this->nxloc];
  if (this->v1_recvW == NULL)  this->v1_recvW = new double[this->nyloc];
  if (this->v1_recvS == NULL)  this->v1_recvS = new double[this->nxloc];
  if (this->v1_sendE == NULL)  this->v1_sendE = new double[this->nyloc];
  if (this->v1_sendN == NULL)  this->v1_sendN = new double[this->nxloc];

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

  // create MPI group from full communicator
  MPI_Group comm_group;
  ierr = MPI_Comm_group(this->comm, &comm_group);
  check_err(ierr, this->comm, "MPI_Comm_group");

  // create nearest-neighbor MPI send/receive groups
  ierr = MPI_Group_incl(comm_group, 1, &(this->nbW), &(this->groupW));
  check_err(ierr, this->comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &(this->nbE), &(this->groupE));
  check_err(ierr, this->comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &(this->nbS), &(this->groupS));
  check_err(ierr, this->comm, "MPI_Group_incl");

  ierr = MPI_Group_incl(comm_group, 1, &(this->nbN), &(this->groupN));
  check_err(ierr, this->comm, "MPI_Group_incl");

  // create RMA window into my receive buffers
  ierr = MPI_Win_create(this->v1_recvW, this->nyloc*sizeof(double), sizeof(double),
                        MPI_INFO_NULL, this->comm, &(this->win_v1W));
  check_err(ierr, this->comm, "MPI_Win_create");

  ierr = MPI_Win_create(this->v1_recvS, this->nxloc*sizeof(double), sizeof(double),
			MPI_INFO_NULL, this->comm, &(this->win_v1S));
  check_err(ierr, this->comm, "MPI_Win_create");

  ierr = MPI_Win_create(this->v2_recvE, this->nyloc*sizeof(double), sizeof(double),
			MPI_INFO_NULL, this->comm, &(this->win_v2E));
  check_err(ierr, this->comm, "MPI_Win_create");

  ierr = MPI_Win_create(this->v3_recvN, this->nxloc*sizeof(double), sizeof(double),
			MPI_INFO_NULL, this->comm, &(this->win_v3N));
  check_err(ierr, this->comm, "MPI_Win_create");

  // free general communicator group
  ierr = MPI_Group_free(&comm_group);
  check_err(ierr, this->comm, "MPI_Group_free");

}


void parallel_decomp::start_v2v3() {

  // expose all east and north recv buffers to others
  int ierr = MPI_Win_post(this->groupE, 0, this->win_v2E);
  check_err(ierr, this->comm, "MPI_Win_post");
  ierr = MPI_Win_post(this->groupN, 0, this->win_v3N);
  check_err(ierr, this->comm, "MPI_Win_post");

  // start RMA access epochs for west and south neighbors' buffers
  ierr = MPI_Win_start(this->groupW, 0, this->win_v2E);
  check_err(ierr, this->comm, "MPI_Win_start");
  ierr = MPI_Win_start(this->groupS, 0, this->win_v3N);
  check_err(ierr, this->comm, "MPI_Win_start");
}


void parallel_decomp::put_v2v3(double *v2, double *v3) {

  // fill send buffers
  for (int j=0; j<this->nyloc; j++)
    this->v2_sendW[j] = v2[idx(0,j,this->nxloc)];
  for (int i=0; i<this->nxloc; i++)
    this->v3_sendS[i] = v3[idx(i,0,this->nxloc)];

  // start MPI_Put calls to put v2_sendW and v3_sendS into neighbor buffers
  int ierr = MPI_Put(this->v2_sendW, this->nyloc, MPI_DOUBLE, this->nbW,
                     0, this->nyloc, MPI_DOUBLE, this->win_v2E);
  check_err(ierr, this->comm, "MPI_Put");

  ierr = MPI_Put(this->v3_sendS, this->nxloc, MPI_DOUBLE, this->nbS,
                 0, this->nxloc, MPI_DOUBLE, this->win_v3N);
  check_err(ierr, this->comm, "MPI_Put");
}


void parallel_decomp::finish_v2v3() {

  // complete the RMA access epochs for neighbor buffers
  int ierr = MPI_Win_complete(this->win_v2E);
  check_err(ierr, this->comm, "MPI_Win_complete");
  ierr = MPI_Win_complete(this->win_v3N);
  check_err(ierr, this->comm, "MPI_Win_complete");

  // wait to ensure that our own v2_recvE and v3_recvN buffers have been filled
  ierr = MPI_Win_wait(this->win_v2E);
  check_err(ierr, this->comm, "MPI_Win_wait");
  ierr = MPI_Win_wait(this->win_v3N);
  check_err(ierr, this->comm, "MPI_Win_wait");
}


void parallel_decomp::start_v1() {

  // expose all west and south recv buffers to others
  int ierr = MPI_Win_post(this->groupW, 0, this->win_v1W);
  check_err(ierr, this->comm, "MPI_Win_post");
  ierr = MPI_Win_post(this->groupS, 0, this->win_v1S);
  check_err(ierr, this->comm, "MPI_Win_post");

  // start RMA access epochs for east and north neighbors' buffers
  ierr = MPI_Win_start(this->groupE, 0, this->win_v1W);
  check_err(ierr, this->comm, "MPI_Win_start");
  ierr = MPI_Win_start(this->groupN, 0, this->win_v1S);
  check_err(ierr, this->comm, "MPI_Win_start");
}


void parallel_decomp::put_v1(double *v1) {

  // fill send buffers
  for (int j=0; j<this->nyloc; j++)
    this->v1_sendE[j] = v1[idx(this->nxloc-1,j,this->nxloc)];
  for (int i=0; i<this->nxloc; i++)
    this->v1_sendN[i] = v1[idx(i,this->nyloc-1,this->nxloc)];

  // start MPI_Put calls to send v1_sendE and v1_sendN
  int ierr = MPI_Put(this->v1_sendE, this->nyloc, MPI_DOUBLE, this->nbE,
                     0, this->nyloc, MPI_DOUBLE, this->win_v1W);
  check_err(ierr, this->comm, "MPI_Put");

  ierr = MPI_Put(this->v1_sendN, this->nxloc, MPI_DOUBLE, this->nbN,
                 0, this->nxloc, MPI_DOUBLE, this->win_v1S);
  check_err(ierr, this->comm, "MPI_Put");
}


void parallel_decomp::finish_v1() {

  // complete the RMA access epochs for neighbor buffers
  int ierr = MPI_Win_complete(this->win_v1W);
  check_err(ierr, this->comm, "MPI_Win_complete");
  ierr = MPI_Win_complete(this->win_v1S);
  check_err(ierr, this->comm, "MPI_Win_complete");

  // wait to ensure that our own v1_recvW and v1_recvS buffers have been filled
  ierr = MPI_Win_wait(this->win_v1W);
  check_err(ierr, this->comm, "MPI_Win_wait");
  ierr = MPI_Win_wait(this->win_v1S);
  check_err(ierr, this->comm, "MPI_Win_wait");
}
