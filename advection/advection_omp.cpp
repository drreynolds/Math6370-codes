/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */


// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "advection.h"
#include "get_time.h"


// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  double stime = get_time();
  /* read problem parameters from input file (should be in this order):
     nx - number of nodes in x-direction
     ny - number of nodes in y-direction
     nt - number of time steps
     tstop - final time (will stop at nt or stop, whichever is 1st)
     c - wave speed
     dtoutput - time frequency for outputting solutions */
  int ierr, nx, ny, nt;
  double tstop, c, dtoutput;
  FILE* FID = fopen("input.txt","r");
  ierr = fscanf(FID," &inputs\n");
  ierr = fscanf(FID,"  nx = %i,\n", &nx);
  ierr = fscanf(FID,"  ny = %i,\n", &ny);
  ierr = fscanf(FID,"  nt = %i,\n", &nt);
  ierr = fscanf(FID,"  tstop = %lf,\n", &tstop);
  ierr = fscanf(FID,"  c = %lf,\n", &c);
  ierr = fscanf(FID,"  dtoutput = %lf,\n", &dtoutput);
  fclose(FID);
  printf("\nRunning wave problem:\n");
  printf("  nx = %i,  ny = %i\n", nx, ny);
  printf("  nt = %i,  tstop = %g\n", nt, tstop);
  printf("  c = %g\n",c);
  printf("  dtoutput = %g\n",dtoutput);
  double ftime = get_time();
  double iotime = ftime-stime;

  // allocate arrays
  stime = get_time();
  double *u  = new double[nx*ny];
  double *v1 = new double[nx*ny];
  double *v2 = new double[nx*ny];
  double *v3 = new double[nx*ny];

  // set grid spacing
  const double dx = 1.0/nx;
  const double dy = 1.0/ny;

  // set initial conditions
  initialize(u, v1, v2, v3, c, dx, dy, nx, ny);
  ftime = get_time();
  double inittime = ftime-stime;

  // set initial time, output initial solution
  double t = 0.0;
  double toutput = t;
  int noutput = 0;
  printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = get_time();
  output(u, toutput, nx, ny, noutput);
  ftime = get_time();
  iotime += ftime-stime;

  // set time step
  const double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start OpenMP parallelism
  double runtime = 0.0;
#pragma omp parallel default(shared)
  {

    // start time stepping
    for (int it=0; it<nt; it++) {

      // start timer
#pragma omp master
      stime = get_time();


      // first update v1 to get to half time step
      // get relevant values for this location
#pragma omp for collapse(2)
      for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {

	  // access relevant components of v2 and v3
          double v2_E = (i == nx-1) ? v2[idx(0,j,nx)] : v2[idx(i+1,j,nx)];
          double v2_W = v2[idx(i,j,nx)];
          double v3_N = (j == ny-1) ? v3[idx(i,0,nx)] : v3_N = v3[idx(i,j+1,nx)];
          double v3_S = v3[idx(i,j,nx)];

	  // update v1
	  v1[idx(i,j,nx)] += c*dt/dx*(v2_E - v2_W)
	    + c*dt/dy*(v3_N - v3_S);

	} // for j
      } // for i (and also end omp for)

      // next update v2 & v3 to get to full time step
      // get relevant values for this location
#pragma omp for collapse(2)
      for (int j=0; j<ny; j++) {
	for (int i=0; i<nx; i++) {

	  // access relevant components of v1
          double v1_W = (i == 0) ? v1[idx(nx-1,j,nx)] : v1[idx(i-1,j,nx)];
          double v1_E = v1[idx(i,j,nx)];
          double v1_S = (j == 0) ? v1[idx(i,ny-1,nx)] : v1_S = v1[idx(i,j-1,nx)];
          double v1_N = v1[idx(i,j,nx)];

	  // update v2 and v3
	  v2[idx(i,j,nx)] += c*dt/dx*(v1_E - v1_W);
	  v3[idx(i,j,nx)] += c*dt/dy*(v1_N - v1_S);

	} // for j
      } // for i (also end omp for)

      // update solution
#pragma omp for collapse(2)
      for (int j=0; j<ny; j++)
	for (int i=0; i<nx; i++) {
	  u[idx(i,j,nx)] += dt*v1[idx(i,j,nx)];
	}  // for i (also end omp for)

      // update time
#pragma omp single
      t = t + dt;
#pragma omp master
      {
	ftime = get_time();
	runtime += ftime-stime;
      }

      // stop simulation if we've reached tstop
      if (t >= tstop)  break;

      // output solution periodically
#pragma omp master
      {
	if (std::abs(t-toutput-dtoutput) <= 1.e-14) {
	  stime = get_time();
	  toutput = t;
	  noutput++;
	  printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
	  output(u, toutput, nx, ny, noutput);
	  ftime = get_time();
	  iotime += ftime-stime;
	}
      }
    } // for it

  }  // end omp parallel region

  // output final solution
  stime = get_time();
  toutput = t;
  noutput++;
  printf("writing output file %i, step = %i, t = %g\n", noutput, nt, t);
  output(u, toutput, nx, ny, noutput);
  ftime = get_time();
  iotime += ftime-stime;

  // output times
  printf(" total initialization time = %.16e\n", inittime);
  printf(" total input/output time   = %.16e\n", iotime);
  printf(" total simulation time     = %.16e\n", runtime);

  // free memory
  delete[] u;
  delete[] v1;
  delete[] v2;
  delete[] v3;

  return 0;
} // end main
