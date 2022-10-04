/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Kokkos_Core.hpp>


// Prototypes
KOKKOS_FUNCTION double maxnorm_dev(double *, int);
KOKKOS_INLINE_FUNCTION double f_dev(double, double);
KOKKOS_INLINE_FUNCTION double fx_dev(double, double);
KOKKOS_INLINE_FUNCTION double fy_dev(double, double);
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

  // set some parameters
  int nx = 400;             // search mesh size
  int ny = 400;             // search mesh size

  // initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

  // set the Kokkos execution space and range policy
#if defined(USE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace;
  std::cout << "Running global minimization with Kokkos using OpenMP backend:\n";
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda       ExecSpace;
  std::cout << "Running global minimization with Kokkos using CUDA backend:\n";
#else
  typedef Kokkos::Serial     ExecSpace;
  std::cout << "Running global minimization with Kokkos using Serial backend:\n";
#endif
  typedef Kokkos::RangePolicy<ExecSpace>    RangePol;

  // start timer
  Kokkos::Timer timer;

  // set subinterval widths
  //   Note: we know the minimum is inside the box [-5,5]x[-5,5]
  const double dx = 10.0/(nx-1);
  const double dy = 10.0/(ny-1);

  // perform steepest descent minimization over all points in the mesh
  typedef Kokkos::MinLoc<double, int>::value_type value_type;
  value_type bestval;
  Kokkos::parallel_reduce( "IC_loop", RangePol(0,nx*ny), KOKKOS_LAMBDA
                           (int i, value_type &mybestval) {

    // set mesh location
    int iy = (i-1)/nx;
    int ix = i - iy*nx;
    double pt[] = {-5.0 + (ix-1)*dx, -5.0 + iy*dy};

    // get current function value
    double fval = f_dev(pt[0],pt[1]);

    // perform a steepest descent minimization at this point
    //    int maxits = 100000000;   // maximum iteration count
    int maxits = 10;   // maximum iteration count
    for (int k=1; k<=maxits; k++) {

      // compute gradient of f at this point
      double df[] = {fx_dev(pt[0],pt[1]), fy_dev(pt[0],pt[1])};

      // set the initial linesearch step size
      double gamma = 1.0/sqrt(df[0]*df[0] + df[1]*df[1]);

      // perform back-tracking line search for gamma
      double tstpt[] = {0.0, 0.0};
      double curval;
      for (int l=1; l<=50; l++) {

        // set test point and calculate function value
        tstpt[0] = pt[0] - gamma*df[0];
        tstpt[1] = pt[1] - gamma*df[1];
        curval = f_dev(tstpt[0],tstpt[1]);

      	// if test point successful, exit; otherwise reduce gamma
      	if (curval < fval)
          break;
        else
          gamma *= 0.5;

      } // end for l

      // check for stagnation/convergence
      double normtest[] = {pt[0] - tstpt[0], pt[1] - tstpt[1]};
      if (maxnorm_dev(normtest,2) < 1.0e-13)  break;

      // update point with current iterate
      pt[0] = tstpt[0];
      pt[1] = tstpt[1];
      fval = curval;

    } // end for k

    // if current value is better than best so far, update best
    if (fval < mybestval.val) {
      mybestval.val = fval;
      mybestval.loc = i;
    }

  }, bestval); // end for i

  // Re-do minimization from best initial guess on host to generate final output
  //   set mesh location
  int iy = (bestval.loc-1)/nx;
  int ix = bestval.loc - iy*nx;
  double pt[] = {-5.0 + (ix-1)*dx, -5.0 + iy*dy};

  //   get current function value
  double fval = f(pt[0],pt[1]);

  //   perform a steepest descent minimization at this point
  int maxits = 100000000;   // maximum iteration count
  for (int k=1; k<=maxits; k++) {

    // compute gradient of f at this point
    double df[] = {fx(pt[0],pt[1]), fy(pt[0],pt[1])};

    // set the initial linesearch step size
    double gamma = 1.0/sqrt(df[0]*df[0] + df[1]*df[1]);

    // perform back-tracking line search for gamma
    double tstpt[] = {0.0, 0.0};
    double curval;
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
    double normtest[] = {pt[0] - tstpt[0], pt[1] - tstpt[1]};
    if (maxnorm(normtest,2) < 1.0e-13)  break;

    // update point with current iterate
    pt[0] = tstpt[0];
    pt[1] = tstpt[1];
    fval = curval;
  } // end for k

  // stop timer
  double runtime = timer.seconds();

  // output computed minimum and corresponding point
  std::cout << "  computed minimum = " << std::setprecision(16) << fval << std::endl;
  std::cout << "             point = (" << pt[0] << ", " << pt[1] << ")" << std::endl;
  std::cout << "           runtime = " << std::setprecision(16) << runtime << std::endl;

  }
  Kokkos::finalize();   // note that the views x, y, and z are automatically deleted here

  return 0;
} // end main


/* Function to compute the max norm of an array,
       || v ||_inf
   where the array v has length n */
KOKKOS_FUNCTION double maxnorm_dev(double *v, int n) {
  double result = 0.0;
  for (int i=0; i<n; i++)
    result = (result > abs(v[i])) ? result : abs(v[i]);
  return result;
}
double maxnorm(double *v, int n) {
  double result = 0.0;
  for (int i=0; i<n; i++)
    result = (result > abs(v[i])) ? result : abs(v[i]);
  return result;
}


// integrand function, f(x,y)
KOKKOS_INLINE_FUNCTION double f_dev(double x, double y) {
  return (exp(sin(50.0*x)) + sin(60.0*exp(y)) + sin(70.0*sin(x)) +
	  sin(sin(80.0*y)) - sin(10.0*(x+y)) + 0.25*(x*x+y*y));
}
inline double f(double x, double y) {
  return (exp(sin(50.0*x)) + sin(60.0*exp(y)) + sin(70.0*sin(x)) +
	  sin(sin(80.0*y)) - sin(10.0*(x+y)) + 0.25*(x*x+y*y));
}


// partial wrt x of integrand function, df/dx
KOKKOS_INLINE_FUNCTION double fx_dev(double x, double y) {
  return (exp(sin(50.0*x))*cos(50.0*x)*50.0 +
          cos(70.0*sin(x))*70.0*cos(x) -
          cos(10.0*(x+y))*10.0 + 0.5*x);
}
inline double fx(double x, double y) {
  return (exp(sin(50.0*x))*cos(50.0*x)*50.0 +
          cos(70.0*sin(x))*70.0*cos(x) -
          cos(10.0*(x+y))*10.0 + 0.5*x);
}


// partial wrt y of integrand function, df/dy
KOKKOS_INLINE_FUNCTION double fy_dev(double x, double y) {
  return (cos(60.0*exp(y))*60.0*exp(y) +
          cos(sin(80.0*y))*cos(80.0*y)*80.0 -
          cos(10.0*(x+y))*10.0 + 0.5*y);
}
inline double fy(double x, double y) {
  return (cos(60.0*exp(y))*60.0*exp(y) +
          cos(sin(80.0*y))*cos(80.0*y)*80.0 -
          cos(10.0*(x+y))*10.0 + 0.5*y);
}
