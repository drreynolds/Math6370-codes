! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


double precision function TIMER()
  !-----------------------------------------------------------------
  ! Function to use the best available timer
  !-----------------------------------------------------------------
  double precision :: dtime
  double precision, external :: omp_get_wtime
  call get_time(TIMER)
  !$TIMER = omp_get_wtime()
end function TIMER


program Global_Min
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes a global minimum to the function
  !        f(x,y) = exp(sin(50x)) + sin(60exp(y)) + sin(70sin(x))
  !               + sin(sin(80y)) - sin(10(x+y)) + (x^2+y^2)/4
  !    We perform a local minimization algorithm (steepest descent) 
  !    using a large number of initial iterates, taken by placing a 
  !    relatively fine discretization over the search space.  We 
  !    note that due to the form of the objective function, we know 
  !    that any global minimum must reside in the box [-5,5]x[-5,5].
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none
  integer, parameter :: nx = 100
  integer, parameter :: ny = 100
  integer, parameter :: maxits = 100000000
  integer :: i, ix, iy, k, l
  double precision :: f, fx, fy, x, y, dx, dy, stime, ftime
  double precision :: pt(2), tstpt(2), df(2), bestpt(2)
  double precision :: curval, fval, gamma, bestval
  double precision, external :: TIMER

  !======= Internals ============
  
  ! set the integrand function
  f(x,y) = exp(sin(5.d1*x)) + sin(6.d1*exp(y)) + sin(7.d1*sin(x)) &
         + sin(sin(8.d1*y)) - sin(1.d1*(x+y)) + (x**2+y**2)/4.d0

  ! set the jacobian functions
  fx(x,y) = exp(sin(5.d1*x))*cos(5.d1*x)*5.d1 &
          + cos(7.d1*sin(x))*7.d1*cos(x)      &
          - cos(1.d1*(x+y))*1.d1 + x/2.d0
  fy(x,y) = cos(6.d1*exp(y))*6.d1*exp(y)      &
          + cos(sin(8.d1*y))*cos(8.d1*y)*8.d1 &
          - cos(1.d1*(x+y))*1.d1 + y/2.d0


  ! start timer
  stime = TIMER()

  ! set subinterval widths
  ! (note: we know that the minimum is inside the box [-5,5]^2)
  dx = 10.d0/(nx-1)
  dy = 10.d0/(ny-1)

  ! perform steepest descent minimization over all points in a mesh
  bestval = 1.d12     ! initialize to very large number

  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(ix,iy,pt,fval,k,df,gamma,l,tstpt,curval)
  do i=1,nx*ny

     iy = (i-1)/nx
     ix = i - iy*nx
     pt(1) = -5.d0 + (ix-1)*dx
     pt(2) = -5.d0 + iy*dy
        
     ! get current function value
     fval = f(pt(1),pt(2))
     
     ! perform a steepest descent minimization at this point
     do k = 1,maxits
        
        ! compute gradient of f at this point
        df(1) = fx(pt(1),pt(2))
        df(2) = fy(pt(1),pt(2))

        ! set the initial linesearch step size
        gamma = 1.d0/sqrt(df(1)**2 + df(2)**2)

        ! perform back-tracking line search for gamma
        do l = 1,50

           ! set test point and calculate function value
           tstpt = pt - gamma*df
           curval = f(tstpt(1),tstpt(2))

           ! if test point successful, exit; otherwise reduce gamma
           if (curval < fval) then
              exit
           else 
              gamma = 0.5d0*gamma
           end if
           
        end do
        
        ! check for stagnation/convergence
        if (maxval(abs((pt-tstpt))) < 1e-13)  exit
        
        ! update point with current iterate
        pt = tstpt
        fval = curval
        
     end do
     
     ! if current value is better than best so far, update best
     !$omp critical
     if (fval < bestval) then
        bestpt = pt
        bestval = fval
        print *,' new best-guess has value ',bestval
     end if
     !$omp end critical

  end do
  !$omp end parallel do

  ! stop timer
  ftime = TIMER()

  ! output computed minimum and corresponding point
  write(*,*) ' computed minimum =',bestval
  write(*,*) '            point =',bestpt(1),bestpt(2)
  write(*,*) '          runtime =',ftime-stime
  

end program Global_Min
!=================================================================
