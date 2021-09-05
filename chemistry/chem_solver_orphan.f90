! Fortran90: ! -*- Mode; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 6370
! 27 March 2011
!=================================================================


subroutine chem_solver_orphan(n,T,u,v,w,lam,eps,maxit,mycount,err)
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
  !    n - integer (input), length of {T, u, v, w} arrays
  !    T - double (input), temperature
  !    u - double (in/out), concentration (in: guess, out: solution)
  !    v - double (in/out), concentration (in: guess, out: solution)
  !    w - double (in/out), concentration (in: guess, out: solution)
  !  lam - double (in), damping parameter (lambda)
  !  eps - double (in), nonlinear solver tolerance (epsilon)
  !  maxit - integer (in), maximum allowed iterations
  !  mycount - integer (out), number of intervals for this thread
  !  err - integer (out), output flag denoting number of failed intervals
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none
  integer, intent(in) :: n, maxit
  integer, intent(out) :: mycount, err
  double precision, intent(in) :: T(n), lam, eps
  double precision, intent(inout) :: u(n), v(n), w(n)
  double precision :: k1, k2, k3, k4, f(3), res
  integer :: i, its

  !======= Internals ============

  ! initialize error flag to success
  err = 0

  ! call solver over n intervals
  mycount = 0
  !$omp do
  do i=1,n,1

     ! compute chemical rate coefficients
     k1 = exp(-5.d0*T(i))
     k2 = atan(5.d0*(T(i)-0.5d0))/3.d0 + 0.5d0
     k3 = 1.d0/cosh(5.d0*(T(i)-0.5d0))
     k4 = tanh(5.d0*(T(i)-0.5d0)**2)
     
     ! compute initial residual function
     f(1) = k2*v(i) + k1*u(i)*(u(i)+v(i)+w(i)-1.d0)
     f(2) = k1*u(i)*(1.d0-u(i)-v(i)-w(i)) - (k2+k3)*v(i) + k4*w(i)
     f(3) = k3*v(i) - k4*w(i)
     res = maxval(abs(f))
     
     ! perform fixed-point iteration
     do its = 0,maxit
        if (res < eps)  exit
        
        ! fixed-point iteration solution update
        u(i) = u(i) + lam*f(1)
        v(i) = v(i) + lam*f(2)
        w(i) = w(i) + lam*f(3)
        
        ! compute residuals
        f(1) = k2*v(i) + k1*u(i)*(u(i)+v(i)+w(i)-1.d0)
        f(2) = k1*u(i)*(1.d0-u(i)-v(i)-w(i)) - (k2+k3)*v(i) + k4*w(i)
        f(3) = k3*v(i) - k4*w(i)
        res = maxval(abs(f))
        
     enddo
     
     if (res < eps) then
!!$        write(*,'(2(A,i8))') '    i =',i,',  its =',its
     else
        write(*,'(2(A,i8),4(A,es9.2))') '    error: i =',i,',  its =',its, &
             ',  res =',res,',  u =',u(i),',  v =',v(i),',  w =',w(i)
        err = err + 1
     endif
     
     ! update thread increment
     mycount = mycount + 1
     
  enddo
  !$omp end do

end subroutine chem_solver_orphan
!=================================================================
