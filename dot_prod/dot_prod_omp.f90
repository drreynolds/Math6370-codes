! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
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


program DotProd_OpenMP
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the dot-product of two vectors.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer, parameter :: MAX_INTEGER_DIGITS = 10
  integer :: n, i, numprocs
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


  ! alloctate the vectors
  stime = TIMER()
  allocate(a(n),b(n))
  ftime = TIMER()
  alloctime = ftime-stime

  ! initialize the vectors
  stime = TIMER()
  do i=1,n
     a(i) = 1.d-3*i/n
     b(i) = 1.d-3*(n-i)/n
  enddo
  ftime = TIMER()
  inittime = ftime-stime

  ! start timer
  stime = TIMER()

  ! start OpenMP parallelism
  sum = 0.d0
  !$omp parallel default(shared)

  ! output parallelism information
  !$ numprocs = omp_get_num_threads()
  !$omp single
  !$ write(*,'(A,i3,A)') ' starting OpenMP with ', numprocs,' processes'
  !$omp end single

  ! perform dot-product 
  !$omp do schedule(static, n/numprocs) reduction(+:sum)
  do i=1,n,1
     sum = sum + a(i)*b(i)
  enddo
  !$omp end do
  !$omp end parallel

  ! stop timer
  ftime = TIMER()
  runtime = ftime-stime

  ! output computed value and runtime
  print *, ' vector length =',n
  print *, '   dot-product =',sum
  print *, '    alloc time =',alloctime
  print *, '     init time =',inittime
  print *, '      run time =',runtime

  ! free vectors
  deallocate(a,b)
  
end program DotProd_OpenMP
!=================================================================
