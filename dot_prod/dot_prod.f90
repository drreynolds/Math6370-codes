! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================


program DotProd
  !-----------------------------------------------------------------
  ! Description:
  !    Computes the dot-product of two vectors.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer, parameter :: MAX_INTEGER_DIGITS = 10
  integer :: n, i
  character(MAX_INTEGER_DIGITS) :: n_string
  double precision, allocatable :: a(:), b(:)
  double precision :: sum, stime, ftime, alloctime, inittime, runtime

  !======= Internals ============

  ! get n from the command line
  call getarg(1, n_string)

  ! ensure that an argument was passed in
  if ( n_string == '' ) then
     stop 'Error: function requires one argument (vector length)'
  endif

  ! convert n_string to integer, and ensure it's positive
  read (n_string, *) n
  if ( n < 1 ) then
     stop 'Error: vector length must be greater than 0'
  endif

  ! allocate the vectors
  call get_time(stime)
  allocate(a(n),b(n))
  call get_time(ftime)
  alloctime = ftime-stime

  ! initialize the vector values
  call get_time(stime)
  do i=1,n
     a(i) = 1.d-3*i/n
     b(i) = 1.d-3*(n-i)/n
  enddo
  call get_time(ftime)
  inittime = ftime-stime

  ! compute dot-product
  call get_time(stime)
  sum = 0.d0
  do i=1,n,1
     sum = sum + a(i)*b(i)
  enddo
  call get_time(ftime)
  runtime = ftime-stime

  ! output computed value and runtime
  print *, ' vector length =',n
  print *, '  dot-product =',sum
  print *, '    alloctime =',alloctime
  print *, '     inittime =',inittime
  print *, '      runtime =',runtime

  ! free vectors
  deallocate(a,b)

end program DotProd
!=================================================================
