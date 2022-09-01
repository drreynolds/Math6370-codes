/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include "get_time.h"

// Prototypes
double maxnorm(double *, int);
double maxval(double *, int);
int maxloc(double *, int);
inline double f(double, double);
inline double fx(double, double);
inline double fy(double, double);


/* Example routine to compute a global minimum to the function
          f(x,y) = exp(sin(50x)) + sin(60exp(y)) + sin(70sin(x))
                 + sin(sin(80y)) - sin(10(x+y)) + (x^2+y^2)/4
   We start with a simple search algorithm that quickly finds 100
   suitable starting points for local minimization algorithms (steepest
   descent).  Each of these starting points are then examined thoroughly
   to find the nearest local minimum. */
int main(int argc, char* argv[]) {

  // set some parameters
  int nx = 1000;            // search mesh size
  int ny = 1000;            // search mesh size
  int npts = 100;           // length of "best" points list
  int maxits = 100000000;   // maximum iteration count

  // start timer
  double stime = get_time();

  // set subinterval widths
  //   Note: we know the minimum is inside the box [-5,5]x[-5,5]
  double dx = 10.0/(nx-1);
  double dy = 10.0/(ny-1);

  printf("initial search over mesh\n");

  // check initial mesh, saving best npts points
  double cutoff = 1.0e12;     // initialize to very large number
  int np = 0;              // no points in queue yet
  double pt[2];
  double curval;
  double *searchpts[2], *searchvals;
  searchpts[0] = new double[npts];
  searchpts[1] = new double[npts];
  searchvals   = new double[npts];
  for (int i=1; i<=nx*ny; i++) {

    // set mesh location
    int iy = (i-1)/nx;
    int ix = i - iy*nx;
    pt[0] = -5.0 + (ix-1)*dx;
    pt[1] = -5.0 + iy*dy;

    // evaluate current point
    curval = f(pt[0],pt[1]);

    // if current value is below cutoff, add to list of points
    if (curval < cutoff) {

      // if list has room, just add this point
      if (np < npts-1) {
        searchpts[0][np] = pt[0];           // add point to list
        searchpts[1][np] = pt[1];
        searchvals[np] = curval;            // add value at point
        np++;
      }
      // if this is the last empty slot in the list, add point
      // and set cutoff
      else if (np == npts-1) {
        searchpts[0][np] = pt[0];           // add point to list
        searchpts[1][np] = pt[1];
        searchvals[np] = curval;            // add value at point
        np++;
        cutoff = maxval(searchvals, npts);  // set cutoff as worst value
      }
      // otherwise replace an inferior entry and update cutoff
      else {
        int idx = maxloc(searchvals, npts); // index of worst pt in list
        searchpts[0][idx] = pt[0];          // replace point to list
        searchpts[1][idx] = pt[1];
        searchvals[idx] = curval;           // replace value
        cutoff = maxval(searchvals, npts);  // update cutoff
      } // end if/else if/else

    } // end if curval

  } // end for i

  printf("performing minimizations over best %i points\n",npts);

  /* We have our list of the best npts test points.  We now do
     local minimization (via steepest descent) around each of
     these to better refine */
  double bestval = 1.0e12;     // initialize to very large number
  double bestpt[2], tstpt[2];
  for (int i=0; i<npts; i++) {

    // extract the current point and its function value
    pt[0] = searchpts[0][i];
    pt[1] = searchpts[1][i];
    double fval = searchvals[i];

    // perform a steepest descent minimization at this point
    for (int k=1; k<=maxits; k++) {

      // compute gradient of f at this point
      double df[2];
      df[0] = fx(pt[0],pt[1]);
      df[1] = fy(pt[0],pt[1]);

      // set the initial linesearch step size
      double gamma = 1.0/sqrt(df[0]*df[0] + df[1]*df[1]);

      // perform back-tracking line search for gamma
      for (int l=1; l<=50; l++) {

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
      double normtest[2];
      normtest[0] = pt[0] - tstpt[0];
      normtest[1] = pt[1] - tstpt[1];
      if (maxnorm(normtest,2) < 1.0e-13)  break;

      // update point with current iterate
      pt[0] = tstpt[0];
      pt[1] = tstpt[1];
      fval = curval;

    } // end for k

    // if current value is better than "best" so far, update best
    if (fval < bestval) {
      bestpt[0] = pt[0];
      bestpt[1] = pt[1];
      bestval = fval;
      printf("  new best-guess has value  %.16e\n",bestval);
    }

  } // end for i

  // stop timer
  double ftime = get_time();
  double runtime = ftime - stime;

  // output computed minimum and corresponding point
  printf("  computed minimum = %.16e\n", bestval);
  printf("             point = (%.11e, %.11e)\n", bestpt[0], bestpt[1]);
  printf("           runtime = %.16e\n", runtime);

  return 0;
} // end main


// Function to find the index of the maximum entry in v (that has length n)
int maxloc(double *v, int n) {
  double mval = v[0];
  int result = 0;
  for (int i=0; i<n; i++)
    if (v[i] > mval) {
      mval = v[i];
      result = i;
    }
  return result;
}


// Function to find the maximum entry in v (that has length n)
double maxval(double *v, int n) {
  double result = v[0];
  for (int i=0; i<n; i++)
    result = (result > v[i]) ? result : v[i];
  return result;
}


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
