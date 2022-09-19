! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
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
  integer :: n, i, numprocs
  character(MAX_INTEGER_DIGITS) :: n_string
  double precision :: a, stime, ftime, zmax, zmaxl
  double precision, allocatable :: x(:), y(:), z(:)
  double precision, external :: omp_get_wtime
  integer, external :: omp_get_num_threads

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
  !$ stime = omp_get_wtime()    ! replace stime with OpenMP version (if available)

  ! start OpenMP parallelism
  a = -3.d0
  zmax = -1.0d-300
  !$omp parallel default(shared) private(i,zmaxl)

  ! ouptut parallelism information (only if OpenMP is enabled)
  !$ numprocs = omp_get_num_threads();
  !$omp single
  !$ print '(A,i2,A)',' starting OpenMP with ',numprocs,' processes'
  !$omp end single

  ! initialize x and y
  !$omp do
  do i=1,n
     x(i) = exp(2.d0*i/n)
     y(i) = 1.d0*(n-1)/n
  end do
  !$omp end do

  ! perform linear combination
  zmaxl = -1.0d300
  !$omp do
  do i=1,n
     z(i) = a*x(i) + y(i)
     x(i) = y(i)/a - z(i)
     y(i) = x(i)*y(i)/n
     zmaxl = max(zmaxl, z(i))
  end do
  !$omp end do

  ! combine maximum values for each thread to global max
  !$omp critical(overall_max)
  zmax = max(zmax, zmaxl)
  !$omp end critical(overall_max)

  ! end parallel region
  !$omp end parallel

  ! output maximum value in z
  write(*,*) '  max(z) =',maxval(z)

  ! stop timer
  call get_time(ftime)
  !$ ftime = omp_get_wtime()

  ! output total time
  write(*,*) ' runtime =',ftime-stime
  
  ! deallocate arrays
  deallocate(x, y, z)


end program axpy
!=================================================================
