/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */


// Inclusions
#include <stdlib.h>
#include <cmath>
#include "advection_mpi.hpp"

// Sets the initial conditions into u, v1, v2, v3
void initialize(double* u, double* v1, double* v2, double* v3, double c,
		double dx, double dy, parallel_decomp& p2d) {


  // allocate temporary arrays
  double *xspan_c = new double[p2d.nxloc];
  double *xspan_h = new double[p2d.nxloc];
  double *yspan_c = new double[p2d.nyloc];
  double *yspan_h = new double[p2d.nyloc];

  // set mesh points
  for (int i=p2d.is; i<=p2d.ie; i++) {
    xspan_c[i-p2d.is] = dx*(0.5 + i);
    xspan_h[i-p2d.is] = dx*i;
  }
  for (int j=p2d.js; j<=p2d.je; j++) {
    yspan_c[j-p2d.js] = dy*(0.5 + j);
    yspan_h[j-p2d.js] = dy*j;
  }

  // set initial condition for solution and derivatives
  for (int j=0; j<p2d.nyloc; j++) {

    // y locations
    double y_c = yspan_c[j];
    double y_h = yspan_h[j];

    for (int i=0; i<p2d.nxloc; i++) {

      // x locations
      double x_c = xspan_c[i];
      double x_h = xspan_h[i];

      /* set initial conditions on u_x, u_y [gaussian blob]
            u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2)) */
      u[ idx(i,j,p2d.nxloc)] = exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v1[idx(i,j,p2d.nxloc)] = 0.0;
      v2[idx(i,j,p2d.nxloc)] = -200.0*c*(x_h-1.0/3.0) *
  	         exp( -100.0*( pow(x_h-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v3[idx(i,j,p2d.nxloc)] = -200.0*c*(y_h-0.5) *
   	         exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_h-0.5,2.0) ) );

    } // for j
  } // for i

  // delete temporary arrays
  delete[] xspan_c;
  delete[] xspan_h;
  delete[] yspan_c;
  delete[] yspan_h;

} // end initialize
