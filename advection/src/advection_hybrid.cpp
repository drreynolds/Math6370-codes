/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   12 November 2021 */

// Inclusions
#include "advection.hpp"
#include "omp.h"

// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // start overall calculation timer
  std::chrono::time_point<std::chrono::system_clock> total_stime =
    std::chrono::system_clock::now();

  /* read problem parameters from input file (should be in this order):
        nx - number of nodes in x-direction
	ny - number of nodes in y-direction
	nt - number of time steps
	tstop - final time (will stop at nt or stop, whichever is 1st)
	c - wave speed
	dtoutput - time frequency for outputting solutions */
  std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();
  int nx, ny, nt;
  double tstop, c, dtoutput;
  FILE* FID = fopen("input.txt","r");
  fscanf(FID," &inputs\n");
  fscanf(FID,"  nx = %i,\n", &nx);
  fscanf(FID,"  ny = %i,\n", &ny);
  fscanf(FID,"  nt = %i,\n", &nt);
  fscanf(FID,"  tstop = %lf,\n", &tstop);
  fscanf(FID,"  c = %lf,\n", &c);
  fscanf(FID,"  dtoutput = %lf,\n", &dtoutput);
  fclose(FID);
  std::cout << "\nRunning wave problem:\n";
  std::cout << "  nx = " << nx << ",  ny = " << ny << std::endl;
  std::cout << "  nt = " << nt << ",  tstop = " << tstop << std::endl;
  std::cout << "  c = " << c << std::endl;
  std::cout << "  dtoutput = " << dtoutput << std::endl;
    std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> iotime = ftime - stime;

  // allocate arrays
  stime = std::chrono::system_clock::now();
  double *u_h, *v1_h, *v2_h, *v3_h, *u_d, *v1_d, *v2_d, *v3_d;
  u_h  = new double[nx*ny];
  v1_h = new double[nx*ny];
  v2_h = new double[nx*ny];
  v3_h = new double[nx*ny];
  cudaMalloc((void**)&u_d,  nx*ny*sizeof(double));
  cudaMalloc((void**)&v1_d, nx*ny*sizeof(double));
  cudaMalloc((void**)&v2_d, nx*ny*sizeof(double));
  cudaMalloc((void**)&v3_d, nx*ny*sizeof(double));

  // create RAJA kernel policy
  using xy_kernel_policy = RAJA::KernelPolicy< RAJA::statement::CudaKernel<
    RAJA::statement::For<1, RAJA::cuda_thread_x_loop,
      RAJA::statement::For<0, RAJA::cuda_thread_y_loop,
        RAJA::statement::Lambda<0> > > > >;

  // set grid spacing
  double dx = 1.0/nx;
  double dy = 1.0/ny;

  // set initial conditions
  initialize(u_h, v1_h, v2_h, v3_h, u_d, v1_d, v2_d, v3_d, c, dx, dy, nx, ny);
  ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> inittime = ftime - stime;

  // set time step
  double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // set initial time, output initial solution
  int output_steps = ceil(dtoutput/dt);
  int noutput = 0;
  std::cout << "writing output file " << noutput << ", step = 0, t = 0" << std::endl;
  stime = std::chrono::system_clock::now();
  output(u_h, u_d, 0.0, nx, ny, noutput);
  ftime = std::chrono::system_clock::now();
  iotime += (ftime - stime);

  // ensure exactly 2 OpenMP threads
  omp_set_num_threads(2);

  // create run timer, and global iteration counter
  int numsteps=0;
  std::chrono::duration<double> runtime;

#pragma omp parallel default(shared)
  {
    // determine thread number
    int myid = omp_get_thread_num();

    // output number of threads in the team
#pragma omp master
    std::cout << "Running with " << omp_get_num_threads() << " threads\n";

    // start time stepping
    for (int it=0; it<nt; it++) {

#pragma omp master
      {
        // start timer
        stime = std::chrono::system_clock::now();

        // update global iteration counter
        numsteps += 1;

        // first update v1 to get to half time step
        RAJA::kernel<xy_kernel_policy>(RAJA::make_tuple(RAJA::RangeSegment(0, ny),
                                                        RAJA::RangeSegment(0, nx)),
                                   [=] RAJA_DEVICE (int j, int i) {

          // access relevant components of v2 and v3
          double v2_E, v2_W, v3_N, v3_S;
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
        RAJA::kernel<xy_kernel_policy>(RAJA::make_tuple(RAJA::RangeSegment(0, ny),
                                                        RAJA::RangeSegment(0, nx)),
                                   [=] RAJA_DEVICE (int j, int i) {

          // access relevant components of v1
          double v1_W, v1_E, v1_S, v1_N;
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

        // update runtime
        ftime = std::chrono::system_clock::now();
        runtime += ftime - stime;
      }

      // set current time
      double t = (it+1)*dt;

      // stop simulation if we've reached tstop
      if (t >= tstop)  break;

      // output solution periodically
      if ( (it+1)%output_steps == 0 ) {

#pragma omp barrier
        if (myid == 0) {
          // copy device data to host
          cudaMemcpy( u_h, u_d, nx*ny*sizeof(double), cudaMemcpyDeviceToHost);
        }
#pragma omp barrier
        if (myid == 1) {
          std::chrono::time_point<std::chrono::system_clock> iostime =  std::chrono::system_clock::now();
          noutput++;
          std::cout << "writing output file " << noutput << ", step = "
                    << it << ", t = " << t << std::endl;        
          output(u_h, u_d, t, nx, ny, noutput);
          std::chrono::time_point<std::chrono::system_clock> ioftime = std::chrono::system_clock::now();
          iotime += ioftime - iostime;
        }
      }

    } // for it
  } // omp parallel

  // output final solution
  stime = std::chrono::system_clock::now();
  noutput++;
  std::cout << "writing output file " << noutput << ", step = "
            << numsteps << ", t = " << numsteps*dt << std::endl;
  output(u_h, u_d, nt*dt, nx, ny, noutput);
  ftime = std::chrono::system_clock::now();
  iotime += ftime - stime;

  // stop overall run timer
  std::chrono::time_point<std::chrono::system_clock> total_ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> total_time = total_ftime - total_stime;

  // output times
  std::cout << " total initialization time = " << std::setprecision(16)
            << inittime.count() << std::endl;
  std::cout << " total input/output time   = " << std::setprecision(16)
            << iotime.count() << std::endl;
  std::cout << " total simulation time     = " << std::setprecision(16)
            << runtime.count() << std::endl;
  std::cout << " total overall time        = " << std::setprecision(16)
            << total_time.count() << std::endl;

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
