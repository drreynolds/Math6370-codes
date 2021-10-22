/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   7 February 2015 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "get_time.h"


// Prototypes
double maxnorm(double *, int);
inline double f(double, double);
inline double fx(double, double);
inline double fy(double, double);


/* Example routine to compute a global minimum to the function
          f(x,y) = exp(sin(50x)) + sin(60exp(y)) + sin(70sin(x))
                 + sin(sin(80y)) - sin(10(x+y)) + (x^2+y^2)/4
   We perform a local minimization algorithm (steepest descent) 
   using a large number of initial iterates, taken by placing a 
   relatively fine discretization over the search space.  We 
   note that due to the form of the objective function, we know 
   that any global minimum must reside in the box [-5,5]x[-5,5]. */
int main(int argc, char* argv[]) {

  // local variables
  int nx, ny, maxits, i, k, l, ix, iy;
  double dx, dy, bestval, fval, gamma, curval, runtime;
  double pt[2], tstpt[2], df[2], normtest[2], bestpt[2];
  double stime, ftime;

  // set some parameters
  nx = 100;             // search mesh size
  ny = 100;             // search mesh size
  maxits = 100000000;   // maximum iteration count

  // start timer
  stime = get_time();

  // set subinterval widths
  //   Note: we know the minimum is inside the box [-5,5]x[-5,5]
  dx = 10.0/(nx-1);
  dy = 10.0/(ny-1);

  // perform steepest descent minimization over all points in the mesh
  bestval = 1.0e12;     // initialize to very large number
  for (i=1; i<=nx*ny; i++) {

    // set mesh location
    iy = (i-1)/nx;
    ix = i - iy*nx;
    pt[0] = -5.0 + (ix-1)*dx;
    pt[1] = -5.0 + iy*dy;
        
    // get current function value
    fval = f(pt[0],pt[1]);

    // perform a steepest descent minimization at this point
    for (k=1; k<=maxits; k++) {
        
      // compute gradient of f at this point
      df[0] = fx(pt[0],pt[1]);
      df[1] = fy(pt[0],pt[1]);

      // set the initial linesearch step size
      gamma = 1.0/sqrt(df[0]*df[0] + df[1]*df[1]);

      // perform back-tracking line search for gamma
      for (l=1; l<=50; l++) {

	      // set test point and calculate function value
	      tstpt[0] = pt[0] - gamma*df[0];
	      tstpt[1] = pt[1] - gamma*df[1];
	      curval = f(tstpt[0],tstpt[1]);

      	// if test point successful, exit; otherwise reduce gamma
      	if (curval < fval) 
	        break;
	      else 
	        gamma *= 0.5;

      } // end for l

      // check for stagnation/convergence
      normtest[0] = pt[0] - tstpt[0];
      normtest[1] = pt[1] - tstpt[1];
      if (maxnorm(normtest,2) < 1.0e-13)  break;
        
      // update point with current iterate
      pt[0] = tstpt[0];
      pt[1] = tstpt[1];
      fval = curval;
        
    } // end for k

    // if current value is better than best so far, update best
    if (fval < bestval) {
      bestpt[0] = pt[0];
      bestpt[1] = pt[1];
      bestval = fval;
      printf("  new best-guess has value  %.16e\n", bestval);
    }
     
  } // end for i

  // stop timer
  ftime = get_time();
  runtime = ftime - stime;

  // output computed minimum and corresponding point
  printf("  computed minimum = %.16e\n",bestval);
  printf("             point = (%.11e, %.11e)\n",bestpt[0],bestpt[1]);
  printf("           runtime = %.16e\n",runtime);

  return 0;
} // end main


/* Function to compute the max norm of an array,
       || v ||_inf 
   where the array v has length n */
double maxnorm(double *v, int n) {
  double result = 0.0;
  for (int i=0; i<n; i++)  
    result = (result > fabs(v[i])) ? result : fabs(v[i]);
  return result;
}


// integrand function, f(x,y)
inline double f(double x, double y) {
  return (exp(sin(50.0*x)) + sin(60.0*exp(y)) + sin(70.0*sin(x)) +
	  sin(sin(80.0*y)) - sin(10.0*(x+y)) + 0.25*(x*x+y*y));
}


// partial wrt x of integrand function, df/dx
inline double fx(double x, double y) {
  return (exp(sin(50.0*x))*cos(50.0*x)*50.0 +
          cos(70.0*sin(x))*70.0*cos(x) -
          cos(10.0*(x+y))*10.0 + 0.5*x);
}


// partial wrt y of integrand function, df/dy
inline double fy(double x, double y) {
  return (cos(60.0*exp(y))*60.0*exp(y) +
          cos(sin(80.0*y))*cos(80.0*y)*80.0 -
          cos(10.0*(x+y))*10.0 + 0.5*y);
}
