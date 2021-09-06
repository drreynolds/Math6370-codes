! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


program DotProd_MPI
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the dot-product of two vectors.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi
!  include "mpif.h"  ! include this if mpi module unavailable

  !======= Declarations =========
  implicit none

  integer, parameter :: MAX_INTEGER_DIGITS = 10
  integer :: n, i, is, ie, numprocs, myid, ierr
  character(MAX_INTEGER_DIGITS) :: n_string
  double precision, allocatable :: a(:), b(:)
  double precision :: mysum, sum, stime, ftime
  double precision :: alloctime, inittime, runtime

  !======= Internals ============
  
  ! intialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

  ! get n from the command line
  call getarg(1, n_string)

  ! ensure that an argument was passed in
  if ( n_string == '' ) then
     write(0,*) 'Error: function requires one argument (vector length)'
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  
  ! convert n_string to integer, and ensure it's positive
  read (n_string, *) n
  if ( n < 1 ) then
     write(0,*) 'Error: vector length must be greater than 0'
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif


  ! allocate the vectors (go ahead and put the whole thing on all procs)
  stime = MPI_Wtime()
  allocate(a(n),b(n))
  ftime = MPI_Wtime()
  alloctime = ftime - stime

  ! initialize the vector values
  stime = MPI_Wtime()
  do i=1,n
     a(i) = 1.d-3*i/n
     b(i) = 1.d-3*(n-i)/n
  enddo
  ftime = MPI_Wtime()
  inittime = ftime - stime

  ! start computation timer
  stime = MPI_Wtime()

  ! determine this processor's interval
  is = floor(1.d0*n/numprocs)*myid + 1
  ie = floor(1.d0*n/numprocs)*(myid+1)
  if (myid == numprocs-1) then
     ie = n
  endif

  ! perform dot-product
  mysum = 0.d0
  do i=is,ie,1
     mysum = mysum + a(i)*b(i)
  enddo

  ! root node collects result
  call MPI_Reduce(mysum, sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                  0, MPI_COMM_WORLD, ierr)

  ! stop timer
  ftime = MPI_Wtime()
  runtime = ftime-stime

  ! output computed value and runtime
  if (myid == 0) then
     print *, ' vector length =',n
     print *, '   dot-product =',sum
     print *, '    alloc time =',alloctime
     print *, '     init time =',inittime
     print *, '      run time =',runtime
  endif

  ! free vectors
  deallocate(a,b)

  ! finalize MPI
  call MPI_Finalize(ierr)  

end program DotProd_MPI
!=================================================================
