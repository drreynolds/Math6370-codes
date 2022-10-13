/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */


/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "advection.h"
#include "get_time.h"



/* Example routine to evolve the first-order 2D wave equations in time. */
int main(int argc, char* argv[]) {

  /* declarations */
  int nx, ny, nt, i, j, noutput, it;
  double tstop, c, dtoutput, *u, *v1, *v2, *v3, dx, dy, t, toutput, dt;
  double v1_N, v1_S, v1_E, v1_W, v2_E, v2_W, v3_N, v3_S, stime, ftime;
  double runtime, iotime, inittime;
  FILE* FID;

  stime = get_time();
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
  printf("  c = %g\n",c);
  printf("  dtoutput = %g\n",dtoutput);
  ftime = get_time();
  iotime = ftime-stime;

  /* allocate arrays */
  stime = get_time();
  u  = malloc(nx*ny * sizeof(double));
  v1 = malloc(nx*ny * sizeof(double));
  v2 = malloc(nx*ny * sizeof(double));
  v3 = malloc(nx*ny * sizeof(double));

  /* set grid spacing */
  dx = 1.0/nx;
  dy = 1.0/ny;

  /* set initial conditions */
  initialize(u, v1, v2, v3, c, dx, dy, nx, ny);
  ftime = get_time();
  inittime = ftime-stime;

  /* set initial time, output initial solution */
  t = toutput = 0.0;
  noutput = 0;
  printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = get_time();
  output(u, t, nx, ny, noutput);
  ftime = get_time();
  iotime += ftime-stime;

  /* set time step */
  dt = fmin(dx/c/50.0, dy/c/50.0);

  /* start OpenMP parallelism */
#pragma omp parallel default(shared) private(v2_E,v2_W,v3_N,v3_S,v1_W,v1_E,v1_S,v1_N,it,i,j)
  {

    /* start time stepping */
#pragma omp single nowait
    runtime = 0.0;
    for (it=0; it<nt; it++) {

      /* start timer */
#pragma omp master
      stime = get_time();

      /* first update v1 to get to half time step
	 get relevant values for this location */
#pragma omp for collapse(2)
      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {

          /* access relevant components of v2 and v3 */
	  if (i == nx-1)   v2_E = v2[idx(0,j,nx)];
	  else             v2_E = v2[idx(i+1,j,nx)];

	  v2_W = v2[idx(i,j,nx)];

	  if (j == ny-1)   v3_N = v3[idx(i,0,nx)];
	  else             v3_N = v3[idx(i,j+1,nx)];

	  v3_S = v3[idx(i,j,nx)];

	  /* update v1 */
	  v1[idx(i,j,nx)] += c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);

	} /* for j */
      } /* for i (also end omp for) */

      /* next update v2 & v3 to get to full time step
	 get relevant values for this location */
#pragma omp for collapse(2)
      for (j=0; j<ny; j++) {
	for (i=0; i<nx; i++) {

          /* access relevant components of v1 */
	  if (i == 0)    v1_W = v1[idx(nx-1,j,nx)];
	  else           v1_W = v1[idx(i-1,j,nx)];

	  v1_E = v1[idx(i,j,nx)];

	  if (j == 0)    v1_S = v1[idx(i,ny-1,nx)];
	  else           v1_S = v1[idx(i,j-1,nx)];

	  v1_N = v1[idx(i,j,nx)];

	  /* update v2 and v3 */
	  v2[idx(i,j,nx)] += c*dt/dx*(v1_E - v1_W);
	  v3[idx(i,j,nx)] += c*dt/dy*(v1_N - v1_S);

	} /* for j */
      } /* for i (also end omp for) */

      /* update solution */
#pragma omp for collapse(2)
      for (j=0; j<ny; j++)
	for (i=0; i<nx; i++) {
	  u[idx(i,j,nx)] += dt*v1[idx(i,j,nx)];
        } /* for i (also end omp for) */

      /* update time */
#pragma omp single
      t = t + dt;
#pragma omp master
      {
        ftime = get_time();
        runtime += ftime-stime;
      }

      /* stop simulation if we've reached tstop */
      if (t >= tstop)  break;

      /* output solution periodically */
#pragma omp master
      {
        if (fabs(t-toutput-dtoutput) <= 1.e-14) {
          stime = get_time();
          toutput = t;
          noutput++;
          printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
          output(u, t, nx, ny, noutput);
          ftime = get_time();
          iotime += ftime-stime;
        }
      }
    } /* for it */

  } /* end omp parallel region */

  /* output final solution */
  stime = get_time();
  toutput = t;
  noutput++;
  printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  output(u, t, nx, ny, noutput);
  ftime = get_time();
  iotime += ftime-stime;

  /* output times */
  printf(" total initialization time = %.16e\n", inittime);
  printf(" total input/output time   = %.16e\n", iotime);
  printf(" total simulation time     = %.16e\n", runtime);

  /* free memory */
  free(u);
  free(v1);
  free(v2);
  free(v3);

  return 0;
} /* end main */
