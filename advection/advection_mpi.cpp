/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   11 May 2017 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "advection_mpi.hpp"
#include "mpi.h"


// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // declarations
  int nx, ny, nt, i, j, it, ierr;
  double tstop, c, dtoutput;
  double v1_N, v1_S, v1_E, v1_W, v2_E, v2_W, v3_N, v3_S;

  // intialize MPI
  parallel_decomp p2d;
  ierr = MPI_Init(&argc, &argv);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Init");

  ierr = MPI_Comm_size(MPI_COMM_WORLD, &(p2d.numprocs));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_size");

  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &(p2d.myid));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Comm_rank");


  double stime = MPI_Wtime();
  if (p2d.myid == 0) {
    /* read problem parameters from input file (should be in this order):
       nx - number of nodes in x-direction
       ny - number of nodes in y-direction
       nt - number of time steps
       tstop - final time (will stop at nt or stop, whichever is 1st)
       c - wave speed
       dtoutput - time frequency for outputting solutions */
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
  }

  // root process broadcasts data to others
  int idata[3];
  double ddata[3];
  if (p2d.myid == 0) {
    idata[0] = nx;
    idata[1] = ny;
    idata[2] = nt;
    ddata[0] = tstop;
    ddata[1] = c;
    ddata[2] = dtoutput;
  }
  ierr = MPI_Bcast(idata, 3, MPI_INT, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");
  ierr = MPI_Bcast(ddata, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Bcast");
  nx = idata[0];
  ny = idata[1];
  nt = idata[2];
  tstop = ddata[0];
  c = ddata[1];
  dtoutput = ddata[2];
  double ftime = MPI_Wtime();
  double iotime = ftime - stime;

  // perform parallel decomposition
  stime = MPI_Wtime();
  p2d.periodic[0] = 1;
  p2d.periodic[1] = 1;
  p2d.pdims[0] = 0;
  p2d.pdims[1] = 0;
  ierr = MPI_Dims_create(p2d.numprocs, 2, p2d.pdims);
  check_err(ierr, MPI_COMM_WORLD, "MPI_Dims_create");

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, p2d.pdims, p2d.periodic, 0, &(p2d.comm));
  check_err(ierr, MPI_COMM_WORLD, "MPI_Cart_create");

  ierr = MPI_Cart_coords(p2d.comm, p2d.myid, 2, p2d.pcoords);
  check_err(ierr, p2d.comm, "MPI_Cart_get");

  int is = p2d.pcoords[0]*nx/p2d.pdims[0];
  int ie = (p2d.pcoords[0]+1)*nx/p2d.pdims[0]-1;
  int js = p2d.pcoords[1]*ny/p2d.pdims[1];
  int je = (p2d.pcoords[1]+1)*ny/p2d.pdims[1]-1;
  
  // allocate arrays
  p2d.nxloc = ie-is+1;
  p2d.nyloc = je-js+1;
  double *u  = new double[p2d.nxloc*p2d.nyloc];
  double *v1 = new double[p2d.nxloc*p2d.nyloc];
  double *v2 = new double[p2d.nxloc*p2d.nyloc];
  double *v3 = new double[p2d.nxloc*p2d.nyloc];
  
  // set grid spacing
  double dx = 1.0/nx;
  double dy = 1.0/ny;

  // set initial conditions
  initialize(u, v1, v2, v3, c, dx, dy, is, ie, js, je);
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
  double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start time stepping 
  double runtime = 0.0;
  double commtime = 0.0;
  for (it=0; it<nt; it++) {

    // communicate v2 (EAST/WEST) & v3 (NORTH/SOUTH)
    stime = MPI_Wtime();
    ierr = p2d.Communication1(v2,v3);
    check_err(ierr, p2d.comm, "Communication1");
    ftime = MPI_Wtime();
    commtime += ftime - stime;
    
    // start timer
    stime = MPI_Wtime();

    // first update v1 to get to half time step
    // get relevant values for this location
    for (j=0; j<p2d.nyloc; j++) {
      for (i=0; i<p2d.nxloc; i++) {

	// access relevant components of v2 and v3
	if (i == p2d.nxloc-1)   v2_E = p2d.v2recvE[j];
	else                    v2_E = v2[idx(i+1,j,p2d.nxloc)];

	v2_W = v2[idx(i,j,p2d.nxloc)];

	if (j == p2d.nyloc-1)   v3_N = p2d.v3recvN[i];
	else                    v3_N = v3[idx(i,j+1,p2d.nxloc)];

	v3_S = v3[idx(i,j,p2d.nxloc)];

	// update v1
	v1[idx(i,j,p2d.nxloc)] += c*dt/dx*(v2_E - v2_W)
                                + c*dt/dy*(v3_N - v3_S);

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
	if (i == 0)    v1_W = p2d.v1recvW[j];
	else           v1_W = v1[idx(i-1,j,p2d.nxloc)];

	v1_E = v1[idx(i,j,p2d.nxloc)];

	if (j == 0)    v1_S = p2d.v1recvS[i];
	else           v1_S = v1[idx(i,j-1,p2d.nxloc)];

	v1_N = v1[idx(i,j,p2d.nxloc)];

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
    printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
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

  // finalize MPI
  ierr = MPI_Finalize();

  return 0;
} // end main


//-------------------------------------------------------


int parallel_decomp::Communication1(double *v2, double *v3) {
  
  // allocate send/recv arrays if NULL
  if (v2recvE == NULL)  v2recvE = new double[nyloc];
  if (v3recvN == NULL)  v3recvN = new double[nxloc];
  if (v2sendW == NULL)  v2sendW = new double[nyloc];
  if (v3sendS == NULL)  v3sendS = new double[nxloc];

  // initialize send/recv buffers
  for (int j=0; j<nyloc; j++)  v2sendW[j] = v2[idx(0,j,nxloc)];
  for (int i=0; i<nxloc; i++)  v3sendS[i] = v3[idx(i,0,nxloc)];
  for (int j=0; j<nyloc; j++)  v2recvE[j] = 0.0;
  for (int i=0; i<nxloc; i++)  v3recvN[i] = 0.0;

  // determine process 'neighbors'
  int ierr;
  if (nbE < 0) {
    int nbcoords[2];
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
  }  // end if (nbE < 0)

  // open send/receive channels
  MPI_Request req[4];
  ierr = MPI_Irecv(v2recvE, nyloc, MPI_DOUBLE, nbE, 2, comm, &(req[0]));
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
  
  // allocate send/recv arrays if NULL
  if (v1recvW == NULL)  v1recvW = new double[nyloc];
  if (v1recvS == NULL)  v1recvS = new double[nxloc];
  if (v1sendE == NULL)  v1sendE = new double[nyloc];
  if (v1sendN == NULL)  v1sendN = new double[nxloc];

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
