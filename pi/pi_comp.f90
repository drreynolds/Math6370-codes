! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================

program ComputePi
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
  integer :: n
  integer :: i
  double precision :: pi, h, x, f, a, stime, ftime
  double precision, parameter :: pi_true = 3.14159265358979323846d0

  !======= Internals ============

  ! set the integrand function
  f(a) = 4.d0 / (1.d0 + a*a)

  ! input the number of intervals
  write(*,*) 'Enter the number of intervals (0 quits):'
  read(*,*) n
  if (n < 1) then
     stop
  endif

  ! start timer
  call get_time(stime)

  ! set subinterval width
  h = 1.d0/n

  ! perform integration over n intervals
  pi = 0.d0
  do i=1,n,1
     x = h*(i - 0.5d0)
     pi = pi + h*f(x)
  enddo

  ! stop timer
  call get_time(ftime)

  ! output computed value and error
  write(*,*) ' computed pi =',pi
  write(*,*) '     true pi =',pi_true
  write(*,*) '       error =',pi_true - pi
  write(*,*) '     runtime =',ftime-stime


end program ComputePi
!=================================================================
