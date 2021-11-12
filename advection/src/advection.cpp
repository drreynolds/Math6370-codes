/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   12 November 2021 */

// Inclusions
#include "advection.hpp"

// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // declarations
  int nx, ny, nt, noutput, it;
  double tstop, c, dtoutput, *u_h, *v1_h, *v2_h, *v3_h, *u_d, *v1_d, *v2_d, *v3_d, dx, dy, t, toutput, dt;
  FILE* FID;

  // start timer
  std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();
    
  /* read problem parameters from input file (should be in this order):
        nx - number of nodes in x-direction
	ny - number of nodes in y-direction
	nt - number of time steps
	tstop - final time (will stop at nt or stop, whichever is 1st)
	c - wave speed
	dtoutput - time frequency for outputting solutions */
  FID = fopen("input.txt","r");
  fscanf(FID," &inputs\n");
  fscanf(FID,"  nx = %i,\n", &nx);
  fscanf(FID,"  ny = %i,\n", &ny);
  fscanf(FID,"  nt = %i,\n", &nt);
  fscanf(FID,"  tstop = %lf,\n", &tstop);
  fscanf(FID,"  c = %lf,\n", &c);
  fscanf(FID,"  dtoutput = %lf,\n", &dtoutput);
  fclose(FID);
  printf("\nRunning wave problem:\n");
  printf("  nx = %i,  ny = %i\n", nx, ny);
  printf("  nt = %i,  tstop = %g\n", nt, tstop);
  printf("  c = %g\n",c);
  printf("  dtoutput = %g\n",dtoutput);
    std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> iotime = ftime - stime;

  // allocate arrays
  stime = std::chrono::system_clock::now();
  u_h  = new double[nx*ny];
  v1_h = new double[nx*ny];
  v2_h = new double[nx*ny];
  v3_h = new double[nx*ny];
  cudaMalloc((void**)&u_d,  nx*ny*sizeof(double));
  cudaMalloc((void**)&v1_d, nx*ny*sizeof(double));
  cudaMalloc((void**)&v2_d, nx*ny*sizeof(double));
  cudaMalloc((void**)&v3_d, nx*ny*sizeof(double));

  // create RAJA policies
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  using xy_kernel_policy = RAJA::KernelPolicy< RAJA::statement::CudaKernel< RAJA::statement::For<1, RAJA::cuda_thread_x_loop, RAJA::statement::For<0, RAJA::cuda_thread_y_loop, RAJA::statement::Lambda<0> > > > >;

  // set grid spacing
  dx = 1.0/nx;
  dy = 1.0/ny;

  // set initial conditions
  initialize(u_h, v1_h, v2_h, v3_h, u_d, v1_d, v2_d, v3_d, c, dx, dy, nx, ny);
  ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> inittime = ftime - stime;

  // set initial time, output initial solution
  t = toutput = 0.0;
  noutput = 0;
  printf("writing output file %i, step = %i, t = %g\n", noutput, 0, t);
  stime = std::chrono::system_clock::now();
  output(u_h, u_d, t, nx, ny, noutput);
  ftime = std::chrono::system_clock::now();
  iotime += (ftime - stime);

  // set time step
  dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start time stepping 
  std::chrono::duration<double> runtime;
  for (it=0; it<nt; it++) {

    // start timer
    stime = std::chrono::system_clock::now();

    // first update v1 to get to half time step
    // get relevant values for this location
    RAJA::kernel<xy_kernel_policy>(RAJA::make_tuple(RAJA::RangeSegment(0, ny),
                                                    RAJA::RangeSegment(0, nx)),
                               [=] RAJA_DEVICE (int j, int i) {

      double v2_E, v2_W, v3_N, v3_S;

      // access relevant components of v2 and v3
      if (i == nx-1)   v2_E = v2_d[idx(0,j,nx)];
	    else             v2_E = v2_d[idx(i+1,j,nx)];

      v2_W = v2_d[idx(i,j,nx)];

      if (j == ny-1)   v3_N = v3_d[idx(i,0,nx)];
	    else             v3_N = v3_d[idx(i,j+1,nx)];

      v3_S = v3_d[idx(i,j,nx)];

      // update v1
      v1_d[idx(i,j,nx)] += c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);
    });

    // next update v2 & v3 to get to full time step
    // get relevant values for this location
    RAJA::kernel<xy_kernel_policy>(RAJA::make_tuple(RAJA::RangeSegment(0, ny),
                                                    RAJA::RangeSegment(0, nx)),
                               [=] RAJA_DEVICE (int j, int i) {

      double v1_W, v1_E, v1_S, v1_N;

	    // access relevant components of v1
      if (i == 0)    v1_W = v1_d[idx(nx-1,j,nx)];
	    else           v1_W = v1_d[idx(i-1,j,nx)];

      v1_E = v1_d[idx(i,j,nx)];

      if (j == 0)    v1_S = v1_d[idx(i,ny-1,nx)];
	    else           v1_S = v1_d[idx(i,j-1,nx)];

      v1_N = v1_d[idx(i,j,nx)];

      // update v2 and v3
	    v2_d[idx(i,j,nx)] += c*dt/dx*(v1_E - v1_W);
	    v3_d[idx(i,j,nx)] += c*dt/dy*(v1_N - v1_S);
    });

    // update solution for plotting
    RAJA::kernel<xy_kernel_policy>(RAJA::make_tuple(RAJA::RangeSegment(0, ny),
                                                    RAJA::RangeSegment(0, nx)),
                               [=] RAJA_DEVICE (int j, int i) {
	    u_d[idx(i,j,nx)] += dt*v1_d[idx(i,j,nx)];
    });

    // update time
    t = t + dt;
    ftime = std::chrono::system_clock::now();
    runtime += ftime - stime;

    // stop simulation if we've reached tstop
    if (t >= tstop)  break;

    // output solution periodically
    if ( t - (toutput + dtoutput) > -1.e-14 ) {
      stime = std::chrono::system_clock::now();
      toutput = t;
      noutput++;
      printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
      output(u_h, u_d, t, nx, ny, noutput);
      ftime = std::chrono::system_clock::now();
      iotime += ftime - stime;
    }

  } // for it
  

  // output final solution
  stime = std::chrono::system_clock::now();
  toutput = t;
  noutput++;
  printf("writing output file %i, step = %i, t = %g\n", noutput, it, t);
  output(u_h, u_d, t, nx, ny, noutput);
  ftime = std::chrono::system_clock::now();
  iotime += ftime - stime;

  // output times
  printf(" total initialization time = %.16e\n", inittime.count());
  printf(" total input/output time   = %.16e\n", iotime.count());
  printf(" total simulation time     = %.16e\n", runtime.count());

  // free memory
  delete[] u_h;
  delete[] v1_h;
  delete[] v2_h;
  delete[] v3_h;
  cudaFree(u_d);
  cudaFree(v1_d);
  cudaFree(v2_d);
  cudaFree(v3_d);

  return 0;
} // end main

