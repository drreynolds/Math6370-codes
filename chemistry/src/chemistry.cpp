/* Daniel R. Reynolds
   SMU Mathematics
   Math 4370/6370
   5 November 2021 */

// Inclusions
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <RAJA/RAJA.hpp>


// Prototypes
RAJA_DEVICE void chem_solver(double, double*, double*, double*,
		 double, double, int, int*, double*);


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

  // allocate temperature and solution arrays
  double *T_h = new double[n];
  double *u_h = new double[n];
  double *v_h = new double[n];
  double *w_h = new double[n];
  double *T_d, *u_d, *v_d, *w_d;
  cudaMalloc((void**)&T_d, n*sizeof(double));
  cudaMalloc((void**)&u_d, n*sizeof(double));
  cudaMalloc((void**)&v_d, n*sizeof(double));
  cudaMalloc((void**)&w_d, n*sizeof(double));

  // set the RAJA policies
  using epolicy = RAJA::cuda_exec<256>;
  using rpolicy = RAJA::cuda_reduce;
  std::cout << "Running chemistry solver with " << n << " intervals, using RAJA with CUDA backend:\n";

  // start timer
  std::chrono::time_point<std::chrono::system_clock> stime =
    std::chrono::system_clock::now();

  // set random temperature field, initial guesses at chemical densities on host
  for (int i=0; i<n; i++)  T_h[i] = random() / (pow(2.0,31.0) - 1.0);
  for (int i=0; i<n; i++)  u_h[i] = 0.35;
  for (int i=0; i<n; i++)  v_h[i] = 0.1;
  for (int i=0; i<n; i++)  w_h[i] = 0.5;

  // transfer T, u, v, w to device
  cudaMemcpy( T_d, T_h, n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy( u_d, u_h, n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy( v_d, v_h, n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy( w_d, w_h, n*sizeof(double), cudaMemcpyHostToDevice);

  // stop initialization phase timer
  std::chrono::time_point<std::chrono::system_clock> ftime =
    std::chrono::system_clock::now();
  std::chrono::duration<double> inittime = ftime - stime;

  // start timer
  stime =
    std::chrono::system_clock::now();

  // call solver over n intervals
  RAJA::ReduceSum<rpolicy, int> errs(0);
  RAJA::forall<epolicy>(RAJA::RangeSegment(0,n), [=] RAJA_DEVICE (int i) {
    int maxit = 1000000;
    double lam = 1.e-2;
    double eps = 1.e-10;
    int its;
    double res;
    chem_solver(T_d[i], &(u_d[i]), &(v_d[i]), &(w_d[i]), lam, eps, maxit, &its, &res);
    if (res >= eps)  errs += 1;
  });

  // copy results back to host
  cudaMemcpy( T_h, T_d, n*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy( u_h, u_d, n*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy( v_h, v_d, n*sizeof(double), cudaMemcpyDeviceToHost);
  cudaMemcpy( w_h, w_d, n*sizeof(double), cudaMemcpyDeviceToHost);

  // stop timer
  ftime = std::chrono::system_clock::now();
  std::chrono::duration<double> runtime = ftime - stime;

  // output solution time
  std::cout << "   inittime = " << std::setprecision(16) << inittime.count() << std::endl;
  std::cout << "    runtime = " << std::setprecision(16) << runtime.count() << std::endl;
  std::cout << "   failures = " << errs.get() << std::endl;

  // delete temperature and solution arrays on host and device
  delete[] T_h;
  delete[] u_h;
  delete[] v_h;
  delete[] w_h;
  cudaFree(T_d);
  cudaFree(u_d);
  cudaFree(v_d);
  cudaFree(w_d);

  return 0;
} // end main



/* Function to compute the equilibrium residuals
          k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
          k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
          k3*v - k4*w                        = fw(u,v,w). */
RAJA_DEVICE void chem_residual(double u, double v, double w, double k1,
		   double k2, double k3, double k4, double f[]) {
  f[0] = k2*v + k1*u*(u+v+w-1.0);
  f[1] = k1*u*(1.0-u-v-w) - (k2+k3)*v + k4*w;
  f[2] = k3*v - k4*w;
}

/* Function to compute the max norm of an array,
       || v ||_inf
   where the array v has length n */
RAJA_DEVICE double maxnorm(double *v, int n) {
  double result=0.0;
  for (int i=0; i<n; i++)
    result = fmax(result, fabs(v[i]));
  return result;
}

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
RAJA_DEVICE void chem_solver(double T, double *u, double *v, double *w, double lam,
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
