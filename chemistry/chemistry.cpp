/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370 / 6370 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <Kokkos_Core.hpp>


// Prototypes
KOKKOS_FUNCTION void chem_solver(double, double*, double*, double*,
                                 double, double, int, int*, double*);
KOKKOS_FUNCTION void chem_residual(double u, double v, double w, double k1,
                                   double k2, double k3, double k4, double f[]);
KOKKOS_FUNCTION double maxnorm(double *v, int n);


/* Example routine to compute the equilibrium chemical densities at
   a number of spatial locations, given a (random) background temperature
   field.  The chemical rate equations and solution strategy are in the
   subroutine chem_solver, which is called at every spatial location. */
int main(int argc, char* argv[]) {

  // ensure that an argument was passed in
  if (argc < 2) {
    std::cerr << "Error: function requires one argument (vector length)\n";
    return 1;
  }

  // set n as the input argument, and ensure it's positive
  int n = atoi(argv[1]);
  if (n < 1) {
    std::cerr << "Error: vector length " << n << " must be greater than 0\n";
    return 1;
  }

  // initialize Kokkos
  Kokkos::initialize( argc, argv );
  {

  // set the Kokkos execution space and range policy
#if defined(USE_OPENMP)
  typedef Kokkos::OpenMP     ExecSpace;
  typedef Kokkos::HostSpace  MemSpace;
  std::cout << "Running chemistry solver with " << n << " intervals, using Kokkos with OpenMP backend:\n";
#elif defined(USE_CUDA)
  typedef Kokkos::Cuda       ExecSpace;
  typedef Kokkos::CudaSpace  MemSpace;
  std::cout << "Running chemistry solver with " << n << " intervals, using Kokkos with CUDA backend:\n";
#elif defined(USE_UVM)
  typedef Kokkos::Cuda          ExecSpace;
  typedef Kokkos::CudaUVMSpace  MemSpace;
  std::cout << "Running chemistry solver with " << n << " intervals, using Kokkos with CUDA backend:\n";
#else
  typedef Kokkos::Serial     ExecSpace;
  typedef Kokkos::HostSpace  MemSpace;
  std::cout << "Running chemistry solver with " << n << " intervals, using Kokkos with Serial backend:\n";
#endif
  typedef Kokkos::RangePolicy<ExecSpace>            RangePol;
  typedef Kokkos::View<double*, MemSpace>           VecViewDev;
  typedef Kokkos::View<double*, Kokkos::HostSpace>  VecViewHost;

  // allocate temperature and solution arrays
  Kokkos::Timer timer;
  VecViewDev  T_d( "T_d", n );
  VecViewDev  u_d( "u_d", n );
  VecViewDev  v_d( "v_d", n );
  VecViewDev  w_d( "w_d", n );
#ifndef USE_UVM  
  VecViewHost T_h( "T_h", n );
  VecViewHost u_h( "u_h", n );
  VecViewHost v_h( "v_h", n );
  VecViewHost w_h( "w_h", n );
#endif
  double alloctime = timer.seconds();

  // set random temperature field, initial guesses at chemical densities on host
  timer.reset();
#ifdef USE_UVM
  for (int i=0; i<n; i++)  T_d(i) = random() / (pow(2.0,31.0) - 1.0);
  for (int i=0; i<n; i++)  u_d(i) = 0.35;
  for (int i=0; i<n; i++)  v_d(i) = 0.1;
  for (int i=0; i<n; i++)  w_d(i) = 0.5;
#else
  for (int i=0; i<n; i++)  T_h(i) = random() / (pow(2.0,31.0) - 1.0);
  for (int i=0; i<n; i++)  u_h(i) = 0.35;
  for (int i=0; i<n; i++)  v_h(i) = 0.1;
  for (int i=0; i<n; i++)  w_h(i) = 0.5;  

  // transfer T, u, v, w to device
  Kokkos::deep_copy( T_d, T_h );
  Kokkos::deep_copy( u_d, u_h );
  Kokkos::deep_copy( v_d, v_h );
  Kokkos::deep_copy( w_d, w_h );
#endif  
  double inittime = timer.seconds();

  // call solver over n intervals
  timer.reset();
  int errs = 0;
  const int maxit = 1000000;
  const double lam = 1.e-2;
  const double eps = 1.e-10;
  Kokkos::parallel_reduce( "solves", RangePol(0,n), KOKKOS_LAMBDA (int i, int &myerrs) {
    int its;
    double res;
    chem_solver(T_d[i], &(u_d[i]), &(v_d[i]), &(w_d[i]), lam, eps, maxit, &its, &res);
    if (res >= eps)  myerrs += 1;
  }, errs);

  // wait for asynchronous calculations to complete
  Kokkos::fence();

#ifndef USE_UVM
  // copy results back to host
  Kokkos::deep_copy( T_h, T_d );
  Kokkos::deep_copy( u_h, u_d );
  Kokkos::deep_copy( v_h, v_d );
  Kokkos::deep_copy( w_h, w_d );
#endif  
  double runtime = timer.seconds();

  // output solution time
  std::cout << "  alloctime = " << std::setprecision(16) << alloctime << std::endl;
  std::cout << "   inittime = " << std::setprecision(16) << inittime << std::endl;
  std::cout << "    runtime = " << std::setprecision(16) << runtime << std::endl;
  std::cout << "   failures = " << errs << std::endl;

  }
  Kokkos::finalize();   // the host/device views T, u, v, and w are automatically deleted here

  return 0;

} // end main



/* Function to compute the equilibrium chemical concentrations, given
   a background temperature field, of a simple reaction network:
          u + x -> v  with rate k1
          u + x <- v  with rate k2
          v -> w      with rate k3
          v <- w      with rate k4,
   where we have assumed that the total concentration (x+u+v+w) = 1,
   and where k1(T), k2(T), k3(T), and k4(T) are the temperature-dependent
   coefficient functions,
          k1(T) = exp(-5*T),
          k2(T) = atan(5*(T-1/2))/3 + 1/2,
          k3(T) = 1/cosh(5*(T-1/2)),
          k4(T) = tanh(5*(T-1/2)^2).
   Using the conservation relation, we write the constrained ODE system
            x = 1 - u - v - w,
          u_t = k2*v - k1*u*x,
          v_t = k1*u*x - (k2+k3)*v + k4*w,
          w_t = k3*v - k4*w.
   Inserting the constraint equation into the rate equations, and
   setting the time derivatives equal to zero (to find equilibrium
   concentrations), we have the system
          0 = k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
          0 = k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
          0 = k3*v - k4*w                        = fw(u,v,w),
   where each of the rate coefficients are frozen at the fixed temperature T.

   To solve this system, we call a simple damped fixed-point iteration:
   given an initial guess X0, compute iterates
          Xn = X{n-1} + lambda*f,
   where 0 < lambda <= 1 is the damping parameter.  We compute these
   iterates Xn until |f(Xn)| < epsilon.

  Arguments:
      T - double   (input), temperature
      u - double*  (in/out), concentration (in: guess, out: solution)
      v - double*  (in/out), concentration (in: guess, out: solution)
      w - double*  (in/out), concentration (in: guess, out: solution)
    lam - double   (in), damping parameter (lambda)
    eps - double   (in), nonlinear solver tolerance (epsilon)
  maxit - integer  (in), maximum allowed iterations
    its - integer* (out), # of iterations required for convergence
    res - double*  (out), final nonlinear residual (max norm)  */
KOKKOS_FUNCTION void chem_solver(double T, double *u, double *v, double *w, double lam,
                                 double eps, int maxit, int *its, double *res) {

  // declarations
  double k1, k2, k3, k4, f[3];
  int i;

  // compute chemical rate coefficients
  k1 = exp(-5.0*T);
  k2 = atan(5.0*(T-0.5))/3.0 + 0.5;
  k3 = 1.0/cosh(5.0*(T-0.5));
  k4 = tanh(5.0*(T-0.5)*(T-0.5));

  // compute initial residual function
  chem_residual(*u, *v, *w, k1, k2, k3, k4, f);
  *res = maxnorm(f, 3);

  // perform fixed-point iteration
  for (i=0; i<=maxit; i++) {
    if (*res < eps)  break;

    // compute fixed-point update
    *u = *u + lam*f[0];
    *v = *v + lam*f[1];
    *w = *w + lam*f[2];

    // compute residuals
    chem_residual(*u, *v, *w, k1, k2, k3, k4, f);
    *res = maxnorm(f, 3);
  }
  *its = i;

} // end chem_solver


/* Function to compute the equilibrium residuals
          k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
          k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
          k3*v - k4*w                        = fw(u,v,w). */
KOKKOS_FUNCTION void chem_residual(double u, double v, double w, double k1,
                                   double k2, double k3, double k4, double f[]) {
  f[0] = k2*v + k1*u*(u+v+w-1.0);
  f[1] = k1*u*(1.0-u-v-w) - (k2+k3)*v + k4*w;
  f[2] = k3*v - k4*w;
}


/* Function to compute the max norm of an array,
       || v ||_inf
   where the array v has length n */
KOKKOS_FUNCTION double maxnorm(double *v, int n) {
  double result=0.0;
  for (int i=0; i<n; i++)
    result = fmax(result, fabs(v[i]));
  return result;
}
