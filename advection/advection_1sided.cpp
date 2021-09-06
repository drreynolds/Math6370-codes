/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "advection_1sided.h"


// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // declarations
  int nx, ny, nxl, nyl, nt, i, is, ie, j, js, je, noutput, it, ierr, nprocs, myid;
  int ibuf[3], pdims[2], periods[2], pcoords[2];
  MPI_Comm comm;
  double dbuf[3];
  double tstop, c, dtoutput, *u, *v1, *v2, *v3, dx, dy, t, toutput, dt;
  double v1_N, v1_S, v1_E, v1_W, v2_E, v2_W, v3_N, v3_S;
  double runtime, iotime, inittime, commtime;
  double *v2_recvE, *v3_recvN, *v1_recvW, *v1_recvS;
  double *v2_sendW, *v3_sendS, *v1_sendE, *v1_sendN;
  FILE* FID;
  double stime, ftime;
  MPI_Win win_v1W, win_v1S, win_v2E, win_v3N;
  MPI_Group groupE, groupW, groupN, groupS;

  // intialize MPI
  ierr = MPI_Init(&argc, &argv);
  if (ierr != MPI_SUCCESS) {
    fprintf(stderr," error in MPI_Init = %i\n",ierr);
    return 1;
  }
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_size");
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_rank");

  stime = MPI_Wtime();
  if (myid == 0) {
    /* read problem parameters from input file (should be in this order):
       nx - number of nodes in x-direction
       ny - number of nodes in y-direction
       nt - number of time steps
       tstop - final time (will stop at nt or stop, whichever is 1st)
       c - wave speed
       dtoutput - time frequency for outputting solutions */
    FID = fopen("input.txt","r");
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
    printf("  nprocs = %i\n", nprocs);
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
  ierr = MPI_Bcast(&ibuf, 5, MPI_INT, 0, MPI_COMM_WORLD);
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

  // set up 2D Cartesian communicator, comm
  periods[0] = 1;
  periods[1] = 1;
  pdims[0] = 0;
  pdims[1] = 0;
  ierr = MPI_Dims_create(nprocs, 2, pdims);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Dims_create");

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, 0, &comm);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Cart_create");

  // get this processor's new id and location from comm
  ierr = MPI_Comm_rank(comm, &myid);
  check_err(ierr, comm, "MPI_Comm_rank");
  ierr = MPI_Cart_coords(comm, myid, 2, pcoords);
  check_err(ierr, comm, "MPI_Cart_coords");

  // determine local extents
  is = ((int) (1.0*nx/pdims[0]))*pcoords[0];
  ie = ((int) (1.0*nx/pdims[0]))*(pcoords[0]+1);
  if (pcoords[0] == pdims[0]-1)  ie = nx;
  js = ((int) (1.0*ny/pdims[1]))*pcoords[1];
  je = ((int) (1.0*ny/pdims[1]))*(pcoords[1]+1);
  if (pcoords[1] == pdims[1]-1)  je = ny;
  nxl = ie-is;
  nyl = je-js;

  // allocate arrays
  u  = new double[nxl*nyl];
  v1 = new double[nxl*nyl];
  v2 = new double[nxl*nyl];
  v3 = new double[nxl*nyl];
  
  // set grid spacing
  dx = 1.0/nx;
  dy = 1.0/ny;

  // set initial conditions
  initialize(u, v1, v2, v3, c, dx, dy, is, js, nxl, nyl);
  ftime = MPI_Wtime();
  inittime = ftime - stime;

  // set initial time, output initial solution
  t = toutput = 0.0;
  noutput = 0;
  if (myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = MPI_Wtime();
  output(u, t, nx, ny, nxl, nyl, noutput, comm);
  ftime = MPI_Wtime();
  iotime += ftime - stime;

  // set time step
  dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // allocate send/recv buffers
  v2_recvE = new double[nyl];
  v3_recvN = new double[nxl];
  v1_recvW = new double[nyl];
  v1_recvS = new double[nxl];
  v2_sendW = new double[nyl];
  v3_sendS = new double[nxl];
  v1_sendE = new double[nyl];
  v1_sendN = new double[nxl];

  // create nearest-neighbor groups and RMA windows into my recv buffers
  stime = MPI_Wtime();
  ierr = setup_windows(nxl, nyl, v1_recvW, v1_recvS, v2_recvE, 
		       v3_recvN, comm, win_v1W, win_v1S, win_v2E, 
		       win_v3N, groupE, groupW, groupN, groupS);
  check_err(ierr, comm, "setup_windows");
  ftime = MPI_Wtime();
  commtime = ftime - stime;

  // start time stepping 
  runtime = 0.0;
  for (it=0; it<nt; it++) {

    // start communication timer
    stime = MPI_Wtime();

    // expose all east and north recv buffers to others
    ierr = MPI_Win_post(groupE, 0, win_v2E);
    check_err(ierr, comm, "MPI_Win_post");
    ierr = MPI_Win_post(groupN, 0, win_v3N);
    check_err(ierr, comm, "MPI_Win_post");

    // start RMA access epochs for west and south neighbors' buffers
    ierr = MPI_Win_start(groupW, 0, win_v2E);
    check_err(ierr, comm, "MPI_Win_start");
    ierr = MPI_Win_start(groupS, 0, win_v3N);
    check_err(ierr, comm, "MPI_Win_start");

    // put v2 and v3 values into neighbor buffers
    ierr = put_v2v3(v2, v3, v2_sendW, v3_sendS, nxl, nyl, 
		    comm, win_v2E, win_v3N);
    check_err(ierr, comm, "put_v2v3");
    
    // complete the RMA access epochs for neighbor buffers
    ierr = MPI_Win_complete(win_v2E);
    check_err(ierr, comm, "MPI_Win_complete");
    ierr = MPI_Win_complete(win_v3N);
    check_err(ierr, comm, "MPI_Win_complete");

    // wait to ensure that our own v2_recvE and v3_recvN buffers have been filled
    ierr = MPI_Win_wait(win_v2E);
    check_err(ierr, comm, "MPI_Win_wait");
    ierr = MPI_Win_wait(win_v3N);
    check_err(ierr, comm, "MPI_Win_wait");

    // update communication timer
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    // start computation timer
    stime = MPI_Wtime();

    // first update v1 to get to half time step
    // get relevant values for this location
    for (j=0; j<nyl; j++) {
      for (i=0; i<nxl; i++) {

	// access relevant components of v2 and v3
	if (i == nxl-1)   v2_E = v2_recvE[j];
	else              v2_E = v2[idx(i+1,j,nxl)];

	v2_W = v2[idx(i,j,nxl)];

	if (j == nyl-1)   v3_N = v3_recvN[i];
	else              v3_N = v3[idx(i,j+1,nxl)];

	v3_S = v3[idx(i,j,nxl)];

	// update v1
	v1[idx(i,j,nxl)] += c*dt/dx*(v2_E - v2_W)
	                  + c*dt/dy*(v3_N - v3_S);

      } // for j
    } // for i

    // update computation timer
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    // start communication timer
    stime = MPI_Wtime();

    // expose all west and south recv buffers to others
    ierr = MPI_Win_post(groupW, 0, win_v1W);
    check_err(ierr, comm, "MPI_Win_post");
    ierr = MPI_Win_post(groupS, 0, win_v1S);
    check_err(ierr, comm, "MPI_Win_post");

    // start RMA access epochs for east and north neighbors' buffers
    ierr = MPI_Win_start(groupE, 0, win_v1W);
    check_err(ierr, comm, "MPI_Win_start");
    ierr = MPI_Win_start(groupN, 0, win_v1S);
    check_err(ierr, comm, "MPI_Win_start");

    // put v1 values into neighbor buffers
    ierr = put_v1(v1, v1_sendE, v1_sendN, nxl, nyl, 
		  comm, win_v1W, win_v1S);
    check_err(ierr, comm, "put_v1");

    // complete the RMA access epochs for neighbor buffers
    ierr = MPI_Win_complete(win_v1W);
    check_err(ierr, comm, "MPI_Win_complete");
    ierr = MPI_Win_complete(win_v1S);
    check_err(ierr, comm, "MPI_Win_complete");

    // wait to ensure that our own v1_recvW and v1_recvS buffers have been filled
    ierr = MPI_Win_wait(win_v1W);
    check_err(ierr, comm, "MPI_Win_wait");
    ierr = MPI_Win_wait(win_v1S);
    check_err(ierr, comm, "MPI_Win_wait");

    // update communication timer
    ftime = MPI_Wtime();
    commtime += ftime - stime;

    // start computation timer
    stime = MPI_Wtime();

    // next update v2 & v3 to get to full time step
    // get relevant values for this location
    for (j=0; j<nyl; j++) {
      for (i=0; i<nxl; i++) {

	// access relevant components of v1
	if (i == 0)    v1_W = v1_recvW[j];
	else           v1_W = v1[idx(i-1,j,nxl)];

	v1_E = v1[idx(i,j,nxl)];

	if (j == 0)    v1_S = v1_recvS[i];
	else           v1_S = v1[idx(i,j-1,nxl)];

	v1_N = v1[idx(i,j,nxl)];

	// update v2 and v3
	v2[idx(i,j,nxl)] += c*dt/dx*(v1_E - v1_W);
	v3[idx(i,j,nxl)] += c*dt/dy*(v1_N - v1_S);
	
      } // for j
    } // for i

    // update solution for plotting
    for (j=0; j<nyl; j++) 
      for (i=0; i<nxl; i++) 
	u[idx(i,j,nxl)] += dt*v1[idx(i,j,nxl)];

    // update time
    t = t + dt;

    // update computation timer
    ftime = MPI_Wtime();
    runtime += ftime - stime;

    // stop simulation if we've reached tstop
    if (t >= tstop)  break;

    // output solution periodically
    if (fabs(t-toutput-dtoutput) <= 1.e-14) {
      stime = MPI_Wtime();
      toutput = t;
      noutput++;
      if (myid == 0)
	printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
      output(u, t, nx, ny, nxl, nyl, noutput, comm);
      ftime = MPI_Wtime();
      iotime += ftime - stime;
    }

  } // for it
  

  // output final solution
  stime = MPI_Wtime();
  toutput = t;
  noutput++;
  if (myid == 0)
    printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
  output(u, t, nx, ny, nxl, nyl, noutput, comm);
  ftime = MPI_Wtime();
  iotime += ftime - stime;

  // output times
  if (myid == 0) {
    printf(" total initialization time = %.16e\n", inittime);
    printf(" total input/output time   = %.16e\n", iotime);
    printf(" total communication time  = %.16e\n", commtime);
    printf(" total simulation time     = %.16e\n", runtime);
  }

  // free RMA windows
  MPI_Win_free(&win_v1W);
  MPI_Win_free(&win_v1S);
  MPI_Win_free(&win_v2E);
  MPI_Win_free(&win_v3N);

  // free MPI groups
  MPI_Group_free(&groupE);
  MPI_Group_free(&groupW);
  MPI_Group_free(&groupN);
  MPI_Group_free(&groupS);

  // free memory
  delete[] u;
  delete[] v1;
  delete[] v2;
  delete[] v3;
  delete[] v3_recvN;
  delete[] v2_recvE;
  delete[] v1_recvS;
  delete[] v1_recvW;

  // finalize MPI
  MPI_Finalize();

  return 0;
} // end main

