/* Daniel R. Reynolds
   SMU Mathematics
   Math 6370
   14 January 2013 */

/* Inclusions */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* Prototypes */
void chem_residual(double u, double v, double w, double k1, 
		   double k2, double k3, double k4, double f[]);
double maxnorm(double *v, int n);


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
int chem_solver_orphan(int n, double *T, double *u, double *v, double *w, 
		       double lam, double eps, int maxit, int *mycount) 
{

  /* declarations*/
  int i, its, err;
  double k1, k2, k3, k4, res, f[3];

  /* initialize error flag to success */
  err = 0;

  /* solve system over n intervals */
  *mycount = 0;
  #pragma omp for
  for (i=0; i<n; i++) {

    /* compute chemical rate coefficients */
    k1 = exp(-5.0*T[i]);
    k2 = atan(5.0*(T[i]-0.5))/3.0 + 0.5;
    k3 = 1.0/cosh(5.0*(T[i]-0.5));
    k4 = tanh(5.0*(T[i]-0.5)*(T[i]-0.5));

    /* compute initial residual function */
    chem_residual(u[i], v[i], w[i], k1, k2, k3, k4, f);
    res = maxnorm(f, 3);

    /* perform fixed-point iteration  */
    for (its=0; its<=maxit; its++) {
      if (res < eps)  break;
   
      /* compute fixed-point update */
      u[i] += lam*f[0];
      v[i] += lam*f[1];
      w[i] += lam*f[2];

      /* compute residuals */
      chem_residual(u[i], v[i], w[i], k1, k2, k3, k4, f);
      res = maxnorm(f, 3);
    }

    if (res < eps) {
      /*       printf("    i = %i,  its = %i\n",i,its);*/
    } else {
      printf("    error: i=%i, its=%i, res=%g, u=%g, v=%g, w=%g\n",
	     i, its, res, u[i], v[i], w[i]);
      err += 1;
    }

    /* update thread increment */
    *mycount += 1;

  }  /* end for  (end omp for)*/

  return err;

} /* end chem_solver_orphan */


/* Function to compute the equilibrium residuals
          k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
          k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
          k3*v - k4*w                        = fw(u,v,w). */
void chem_residual(double u, double v, double w, double k1, 
		   double k2, double k3, double k4, double f[])
{
  f[0] = k2*v + k1*u*(u+v+w-1.0);
  f[1] = k1*u*(1.0-u-v-w) - (k2+k3)*v + k4*w;
  f[2] = k3*v - k4*w;
}


/* Function to compute the max norm of an array,
       || v ||_inf 
   where the array v has length n */
double maxnorm(double *v, int n)
{
  double result=0.0;
  int i;
  for (i=0; i<n; i++)  
    result = fmax(result, fabs(v[i]));
  return result;
}
