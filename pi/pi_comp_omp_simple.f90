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


program ComputePi_OpenMP
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes pi through numerical integration via
  !        pi = 4*int_0^1 1/(1+x^2) dx
  !    We use a simple midpoint rule for integration, over 
  !    subintervals of fixed size 1/n, where n is a user-input 
  !    parameter.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer :: n, i
  double precision :: pi, h, x, f, a, stime, ftime
  double precision, parameter :: pi_true = 3.14159265358979323846d0
  double precision, external :: TIMER

  !======= Internals ============
  
  ! set the integrand function
  f(a) = 4.d0 / (1.d0 + a*a)

  ! input the number of intervals
  print *, 'Enter the number of intervals (0 quits):'
  read(*,*) n
  if (n < 1) then
     stop
  endif

  ! start timer
  stime = TIMER()

  ! set subinterval width
  h = 1.d0/n

  ! perform integration over n intervals, parallelizing the loop
  pi = 0.d0
  !$omp parallel do private(x) reduction(+:pi)
  do i=1,n,1
     x = h*(i - 0.5d0)
     pi = pi + h*f(x)
  enddo
  !$omp end parallel do

  ! stop timer
  ftime = TIMER()

  ! output computed value and error
  print *, ' computed pi =', pi
  print *, '     true pi =', pi_true
  print *, '       error =', pi_true-pi
  print *, '     runtime =', ftime-stime
  

end program ComputePi_OpenMP
!=================================================================
