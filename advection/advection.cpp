/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include "advection.hpp"

// Example routine to evolve the first-order 2D wave equations in time
int main(int argc, char* argv[]) {

  // initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

  // start overall calculation timer
  Kokkos::Timer total_timer;

  /* read problem parameters from input file (should be in this order):
        nx - number of nodes in x-direction
	ny - number of nodes in y-direction
	nt - number of time steps
	tstop - final time (will stop at nt or stop, whichever is 1st)
	c - wave speed
	dtoutput - time frequency for outputting solutions */
  Kokkos::Timer timer;
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
#if defined(USE_OPENMP)
  std::cout << "\nRunning wave problem using Kokkos with OpenMP backend:\n";
#elif defined(USE_CUDA)
  std::cout << "\nRunning wave problem using Kokkos with CUDA backend:\n";
#else
  std::cout << "\nRunning wave problem using Kokkos with Serial backend:\n";
#endif
  std::cout << "  nx = " << nx << ",  ny = " << ny << std::endl;
  std::cout << "  nt = " << nt << ",  tstop = " << tstop << std::endl;
  std::cout << "  c = " << c << std::endl;
  std::cout << "  dtoutput = " << dtoutput << std::endl;
  double iotime = timer.seconds();

  // allocate arrays
  timer.reset();
  Vec2D  u_d(  "u_d",  nx, ny );
  Vec2D  v1_d( "v1_d", nx, ny );
  Vec2D  v2_d( "v2_d", nx, ny );
  Vec2D  v3_d( "v3_d", nx, ny );
  Vec2DHost u_h  = Kokkos::create_mirror_view(u_d);
  Vec2DHost v1_h = Kokkos::create_mirror_view(v1_d);
  Vec2DHost v2_h = Kokkos::create_mirror_view(v2_d);
  Vec2DHost v3_h = Kokkos::create_mirror_view(v3_d);

  // set grid spacing
  const double dx = 1.0/nx;
  const double dy = 1.0/ny;

  // set initial conditions
  initialize(u_h, v1_h, v2_h, v3_h, u_d, v1_d, v2_d, v3_d, c, dx, dy, nx, ny);
  double inittime = timer.seconds();

  // set initial time, output initial solution
  double t = 0.0;
  double toutput = 0.0;
  int noutput = 0;
  std::cout << "writing output file " << noutput << ", step = 0, t = " << t << std::endl;
  timer.reset();
  output(u_h, u_d, t, nx, ny, noutput);
  iotime += timer.seconds();

  // set time step
  const double dt = (dx < dy) ? dx/c/50.0 : dy/c/50.0;

  // start time stepping
  int it;
  double runtime;
  for (it=0; it<nt; it++) {

    // start timer
    timer.reset();

    // first update v1 to get to half time step
    Kokkos::parallel_for( "update_v1", DevRange2D({0,0},{nx,ny}), KOKKOS_LAMBDA (int i, int j) {

      // access relevant components of v2 and v3
      double v2_E = (i == nx-1) ? v2_d(0,j) : v2_d(i+1,j);
      double v2_W = v2_d(i,j);
      double v3_N = (j == ny-1) ? v3_d(i,0) : v3_d(i,j+1);
      double v3_S = v3_d(i,j);

      // update v1
      v1_d(i,j) += c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S);
    });
    Kokkos::fence();

    // next update v2 & v3 to get to full time step
    Kokkos::parallel_for( "update_v2_v3", DevRange2D({0,0},{nx,ny}), KOKKOS_LAMBDA (int i, int j) {

      // access relevant components of v1
      double v1_W = (i == 0) ? v1_d(nx-1,j) : v1_d(i-1,j);
      double v1_E = v1_d(i,j);
      double v1_S = (j == 0) ? v1_d(i,ny-1) : v1_d(i,j-1);
      double v1_N = v1_d(i,j);

      // update v2 and v3
      v2_d(i,j) += c*dt/dx*(v1_E - v1_W);
      v3_d(i,j) += c*dt/dy*(v1_N - v1_S);
    });
    Kokkos::fence();

    // update solution for plotting
    Kokkos::parallel_for( "update_u", DevRange2D({0,0},{nx,ny}), KOKKOS_LAMBDA (int i, int j) {
      u_d(i,j) += dt*v1_d(i,j);
    });
    Kokkos::fence();

    // update time
    t += dt;
    runtime += timer.seconds();

    // stop simulation if we've reached tstop
    if (t >= tstop)  break;

    // output solution periodically
    if ( t - (toutput + dtoutput) > -1.e-14 ) {
      timer.reset();
      toutput = t;
      noutput++;
      std::cout << "writing output file " << noutput << ", step = "
                << it << ", t = " << t << std::endl;
      output(u_h, u_d, t, nx, ny, noutput);
      iotime += timer.seconds();
    }

  } // for it

  // output final solution
  timer.reset();
  toutput = t;
  noutput++;
  std::cout << "writing output file " << noutput << ", step = "
            << it << ", t = " << t << std::endl;
  output(u_h, u_d, t, nx, ny, noutput);
  iotime += timer.seconds();

  // stop overall run timer
  double total_time = total_timer.seconds();

  // output times
  std::cout << " total initialization time = " << std::setprecision(16)
            << inittime << std::endl;
  std::cout << " total input/output time   = " << std::setprecision(16)
            << iotime << std::endl;
  std::cout << " total simulation time     = " << std::setprecision(16)
            << runtime << std::endl;
  std::cout << " total overall time        = " << std::setprecision(16)
            << total_time << std::endl;

  }
  Kokkos::finalize();   // the host/device views for {u,v1,v2,v3} are automatically deleted here

  return 0;
} // end main
