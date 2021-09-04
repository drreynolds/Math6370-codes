/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */


// Inclusions 
#include <stdlib.h>     // new, delete
#include <math.h>       // exp(), pow()
#include "advection.h"  // idx(), prototypes

// Sets the initial conditions into u, v1, v2, v3
void initialize(double* u, double* v1, double* v2, double* v3, 
		double c, double dx, double dy, int nx, int ny) {

  // declarations
  double *xspan_c, *xspan_h, *yspan_c, *yspan_h, y_c, y_h, x_c, x_h;
  int i, j;

  // allocate temporary arrays
  xspan_c = new double[nx];
  xspan_h = new double[nx];
  yspan_c = new double[ny];
  yspan_h = new double[ny];
  
  // set mesh points
  for (i=0; i<nx; i++) {
    xspan_c[i] = dx*(0.5 + i);
    xspan_h[i] = dx*i;
  }
  for (j=0; j<ny; j++) {
    yspan_c[j] = dy*(0.5 + j);
    yspan_h[j] = dy*j;
  }

  // set initial condition for solution and derivatives
  for (j=0; j<ny; j++) {

    // y locations
    y_c = yspan_c[j];
    y_h = yspan_h[j];

    for (i=0; i<nx; i++) {

      // x locations
      x_c = xspan_c[i];
      x_h = xspan_h[i];

      /* set initial conditions on u_x, u_y [gaussian blob]
            u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2)) */
      u[ idx(i,j,nx)] = exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v1[idx(i,j,nx)] = 0.0;
      v2[idx(i,j,nx)] = -200.0*c*(x_h-1.0/3.0) *
  	         exp( -100.0*( pow(x_h-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v3[idx(i,j,nx)] = -200.0*c*(y_h-0.5) * 
   	         exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_h-0.5,2.0) ) );

    } // for j
  } // for i

  // delete temporary arrays
  delete[] xspan_c;
  delete[] xspan_h;
  delete[] yspan_c;
  delete[] yspan_h;

} // end initialize

