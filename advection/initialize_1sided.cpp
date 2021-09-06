/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   30 April 2015 */


// Inclusions 
#include <stdlib.h>     // new, delete
#include <math.h>       // exp(), pow()
#include "advection_1sided.h"

// Sets the initial conditions into u, v1, v2, v3
void initialize(double* u, double* v1, double* v2, 
		double* v3, double c, double dx, double dy, 
		int is, int js, int nxl, int nyl) {

  // declarations
  double y_c, y_h, x_c, x_h;
  double *xspan_c, *xspan_h, *yspan_c, *yspan_h;
  int i, j;

  // allocate temporary arrays
  xspan_c = new double[nxl];
  xspan_h = new double[nxl];
  yspan_c = new double[nyl];
  yspan_h = new double[nyl];
  
  // set mesh points
  for (i=is; i<is+nxl; i++) {
    xspan_c[i-is] = dx*(0.5 + i);
    xspan_h[i-is] = dx*i;
  }
  for (j=0; j<nyl; j++) {
    yspan_c[j] = dy*(0.5 + j + js);
    yspan_h[j] = dy*(j + js);
  }

  // set initial condition for solution and derivatives
  for (j=0; j<nyl; j++) {

    // y locations
    y_c = yspan_c[j];
    y_h = yspan_h[j];

    for (i=0; i<nxl; i++) {

      // x locations
      x_c = xspan_c[i];
      x_h = xspan_h[i];

      /* set initial conditions on u_x, u_y [gaussian blob]
            u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2)) */
      u[ idx(i,j,nxl)] = exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v1[idx(i,j,nxl)] = 0.0;
      v2[idx(i,j,nxl)] = -200.0*c*(x_h-1.0/3.0) *
  	         exp( -100.0*( pow(x_h-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v3[idx(i,j,nxl)] = -200.0*c*(y_h-0.5) * 
   	         exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_h-0.5,2.0) ) );

    } // for j
  } // for i

  // delete temporary arrays
  delete[] xspan_c;
  delete[] xspan_h;
  delete[] yspan_c;
  delete[] yspan_h;

} // end initialize

