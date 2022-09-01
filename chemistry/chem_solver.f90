! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================

subroutine chem_solver(T,u,v,w,lam,eps,maxit,its,res)
  !-----------------------------------------------------------------
  ! Description:
  !    Computes the equilibrium chemical concentrations, given a
  !    background temperature field, of a simple reaction network:
  !        u + x -> v  with rate k1
  !        u + x <- v  with rate k2
  !        v -> w      with rate k3
  !        v <- w      with rate k4,
  !    where we have assumed that the total concentration
  !    (x+u+v+w) = 1, and where k1(T), k2(T), k3(T), and k4(T) are
  !    the temperature-dependent coefficient functions,
  !        k1(T) = exp(-5*T),
  !        k2(T) = atan(5*(T-1/2))/3 + 1/2,
  !        k3(T) = 1/cosh(5*(T-1/2)),
  !        k4(T) = tanh(5*(T-1/2)^2).
  !    Using the conservation relation, we write the constrained
  !    ODE system
  !          x = 1 - u - v - w,
  !        u_t = k2*v - k1*u*x,
  !        v_t = k1*u*x - (k2+k3)*v + k4*w,
  !        w_t = k3*v - k4*w.
  !    Inserting the constraint equation into the rate equations,
  !    and setting the time derivatives equal to zero (to find
  !    equilibrium concentrations), we have the system
  !        0 = k2*v + k1*u*(u+v+w-1)              = fu(u,v,w),
  !        0 = k1*u*(1-u-v-w) - (k2+k3)*v + k4*w  = fv(u,v,w),
  !        0 = k3*v - k4*w                        = fw(u,v,w),
  !    where each of the rate coefficients are frozen at the fixed
  !    temperature T.
  !
  !    To solve this system, we call a simple damped fixed-point
  !    iteration: given an initial guess X0, compute iterates
  !        Xn = X{n-1} + lambda*f,
  !    where 0 < lambda <= 1 is the damping parameter.  We compute
  !    these iterates Xn until |f(Xn)| < epsilon.
  !
  ! Arguments:
  !    T - double (input), temperature
  !    u - double (in/out), concentration (in: guess, out: solution)
  !    v - double (in/out), concentration (in: guess, out: solution)
  !    w - double (in/out), concentration (in: guess, out: solution)
  !  lam - double (in), damping parameter (lambda)
  !  eps - double (in), nonlinear solver tolerance (epsilon)
  !  maxit - integer (in), maximum allowed iterations
  !  its - integer (out), # of iterations required for convergence
  !  res - double (out), final nonlinear residual (max norm)
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none
  integer, intent(in) :: maxit
  integer, intent(out) :: its
  double precision, intent(in) :: T, lam, eps
  double precision, intent(out) :: res
  double precision, intent(inout) :: u, v, w

  double precision :: k1, k2, k3, k4, f(3)

  !======= Internals ============

  ! compute chemical rate coefficients
  k1 = exp(-5.d0*T)
  k2 = atan(5.d0*(T-0.5d0))/3.d0 + 0.5d0
  k3 = 1.d0/cosh(5.d0*(T-0.5d0))
  k4 = tanh(5.d0*(T-0.5d0)**2)

  ! compute initial residual function
  f = chem_residual()
  res = maxval(abs(f))

  ! loop
  do its = 0,maxit

     if (res < eps)  exit

     ! fixed-point iteration solution update
     u = u + lam*f(1)
     v = v + lam*f(2)
     w = w + lam*f(3)

     ! compute residuals
     f = chem_residual()
     res = maxval(abs(f))

  enddo

contains

  function chem_residual()
    double precision :: chem_residual(3)
    chem_residual(1) = k2*v + k1*u*(u+v+w-1.d0)
    chem_residual(2) = k1*u*(1.d0-u-v-w) - (k2+k3)*v + k4*w
    chem_residual(3) = k3*v - k4*w
  end function chem_residual

end subroutine chem_solver
!=================================================================
