! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 6370
! 22 February 2009
!=================================================================


program axpy
  !-----------------------------------------------------------------
  ! Description: 
  !    Performs some simple vector linear combinations.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer, parameter :: MAX_INTEGER_DIGITS = 10
  integer :: n, i
  character(MAX_INTEGER_DIGITS) :: n_string
  double precision :: a, stime, ftime
  double precision, allocatable :: x(:), y(:), z(:)

  !======= Internals ============
  
  ! get n from the command line
  call getarg(1, n_string)

  ! ensure that an argument was passed in
  if ( n_string == '' ) then
     stop 'Error: function requires one argument (vector length)'
  endif
  
  ! convert n_string to integer
  read (n_string, *) n

  ! allocate arrays
  allocate(x(n), y(n), z(n))

  ! start timer
  call get_time(stime)

  ! initialize a, x and y
  a = -3.d0
  do i=1,n
     x(i) = exp(2.d0*i/n)
     y(i) = 1.d0*(n-1)/n
  end do

  ! perform linear combination
  do i=1,n
     z(i) = a*x(i) + y(i)
  end do
  
  do i=1,n
     x(i) = y(i)/a - z(i)
  end do
  
  do i=1,n
     y(i) = x(i)*y(i)/n
  end do

  ! output maximum value in z
  write(*,*) '  max(z) =',maxval(z)

  ! stop timer
  call get_time(ftime)

  ! output total time
  write(*,*) ' runtime =',ftime-stime
  
  ! deallocate arrays
  deallocate(x, y, z)


end program axpy
!=================================================================
