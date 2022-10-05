/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include "advection.hpp"

// Sets the initial conditions into u, v1, v2, v3
void initialize(Vec2DHost u_h, Vec2DHost v1_h, Vec2DHost v2_h, Vec2DHost v3_h,
                Vec2D u_d, Vec2D v1_d, Vec2D v2_d, Vec2D v3_d,
                double c, double dx, double dy, int nx, int ny) {

  // set temporary mesh points (in serial, on host)
  Kokkos::View<double*, Kokkos::HostSpace> xspan_c( "x_c", nx );
  Kokkos::View<double*, Kokkos::HostSpace> xspan_h( "x_h", nx );
  for (int i=0; i<nx; i++) {
    xspan_c(i) = dx*(0.5 + i);
    xspan_h(i) = dx*i;
  }
  Kokkos::View<double*, Kokkos::HostSpace> yspan_c( "y_c", ny );
  Kokkos::View<double*, Kokkos::HostSpace> yspan_h( "y_c", ny );
  for (int j=0; j<ny; j++) {
    yspan_c(j) = dy*(0.5 + j);
    yspan_h(j) = dy*j;
  }

  // set initial condition for solution and derivatives
  Kokkos::parallel_for( "set_ICs", HostRange2D({0,0},{nx,ny}), KOKKOS_LAMBDA (int i, int j) {

    // (x,y) locations
    double y_c = yspan_c(j);
    double y_h = yspan_h(j);
    double x_c = xspan_c(i);
    double x_h = xspan_h(i);

    /* set initial conditions on u_x, u_y [gaussian blob]
          u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
          c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
          c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2)) */
    u_h(i,j)  = exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
    v1_h(i,j) = 0.0;
    v2_h(i,j) = -200.0*c*(x_h-1.0/3.0) *
      exp( -100.0*( pow(x_h-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
    v3_h(i,j) = -200.0*c*(y_h-0.5) *
      exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_h-0.5,2.0) ) );

  });
  Kokkos::fence();

  // copy initial states to device
  Kokkos::deep_copy( u_d,   u_h  );
  Kokkos::deep_copy( v1_d,  v1_h );
  Kokkos::deep_copy( v2_d,  v2_h );
  Kokkos::deep_copy( v3_d,  v3_h );

  // locally-declared host views will be automatically deleted upon function return

} // end initialize
