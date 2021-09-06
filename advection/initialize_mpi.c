/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   11 May 2017 */

/* Inclusions */
#include <stdlib.h>
#include <math.h>
#include "advection_mpi.h"
#include "mpi.h"


/* Sets the initial conditions into u, v1, v2, v3. */
void initialize(double* u, double* v1, double* v2, double* v3, double c,
                double dx, double dy, int is, int ie, int js, int je) {

  /* declarations */
  double *xspan_c, *xspan_h, *yspan_c, *yspan_h, y_c, y_h, x_c, x_h;
  int i, j, nxloc=ie-is+1, nyloc=je-js+1;

  /* allocate temporary arrays */
  xspan_c = (double *) malloc(nxloc * sizeof(double));
  xspan_h = (double *) malloc(nxloc * sizeof(double));
  yspan_c = (double *) malloc(nyloc * sizeof(double));
  yspan_h = (double *) malloc(nyloc * sizeof(double));
  
  /* set local coordinates */
  for (i=is; i<=ie; i++) {
    xspan_c[i-is] = dx*(0.5 + i);
    xspan_h[i-is] = dx*i;
  }
  for (j=js; j<=je; j++) {
    yspan_c[j-js] = dy*(0.5 + j);
    yspan_h[j-js] = dy*j;
  }

  /* set initial condition for solution and derivatives */
  for (j=0; j<nyloc; j++) {
    
    /* y locations */
    y_c = yspan_c[j];
    y_h = yspan_h[j];

    for (i=0; i<nxloc; i++) {

      /* x locations */
      x_c = xspan_c[i];
      x_h = xspan_h[i];

      /* set initial conditions on u_x, u_y [gaussian blob]
            u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
            c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2)) */
      u[ idx(i,j,nxloc)]  = exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v1[idx(i,j,nxloc)] = 0.0;
      v2[idx(i,j,nxloc)] = -200.0*c*(x_h-1.0/3.0) *
  	         exp( -100.0*( pow(x_h-1.0/3.0,2.0) + pow(y_c-0.5,2.0) ) );
      v3[idx(i,j,nxloc)] = -200.0*c*(y_h-0.5) * 
   	         exp( -100.0*( pow(x_c-1.0/3.0,2.0) + pow(y_h-0.5,2.0) ) );

    } /* for j */
  } /* for i */

  /* free temporary arrays */
  free(xspan_c);
  free(xspan_h);
  free(yspan_c);
  free(yspan_h);

} /* end initialize_mpi */

