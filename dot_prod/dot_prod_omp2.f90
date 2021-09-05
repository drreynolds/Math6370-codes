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
  double precision, external :: TIMER
  integer, external :: omp_get_num_threads

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

  ! start parallel region
  !$omp parallel default(none) private(i) &
  !$omp shared(a,b,n,stime,ftime,sum,alloctime,inittime,runtime) 

  ! allocate the vectors
  !$omp single
  stime = TIMER()
  allocate(a(n),b(n))
  ftime = TIMER()
  alloctime = ftime-stime
  !$omp end single

  ! initialize the vectors
  !$omp master
  stime = TIMER()
  !$omp end master
  !$omp do
  do i=1,n
     a(i) = 1.d-3*i/n
     b(i) = 1.d-3*(n-i)/n
  enddo
  !$omp end do
  !$omp master
  ftime = TIMER()
  inittime = ftime-stime
  !$omp end master

  ! start timer
  !$omp master
  stime = TIMER()
  !$omp end master

  ! compute dot-product
  !$omp single
  sum = 0.d0
  !$omp end single
  !$omp do reduction(+:sum)
  do i=1,n,1
     sum = sum + a(i)*b(i)
  enddo
  !$omp end do

  ! stop timer
  !$omp master
  ftime = TIMER()
  runtime = ftime-stime
  !$omp end master

  ! output computed value and runtime
  !$omp master
  print *, ' vector length =',n
  !$ print *, ' num threads =', omp_get_num_threads()
  print *, '   dot-product =',sum
  print *, '    alloc time =',alloctime
  print *, '     init time =',inittime
  print *, '      run time =',runtime
  !$omp end master

  ! free vectors
  !$omp single
  deallocate(a,b)
  !$omp end single

  !$omp end parallel

end program DotProd
!=================================================================
