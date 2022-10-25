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
  use mpi
!  include "mpif.h"  ! include this if mpi module unavailable

  !======= Declarations =========
  implicit none

  integer :: n, i, is, ie, numprocs, myid, ierr
  double precision :: mypi, pi, h, x, f, a, stime, ftime
  double precision, parameter :: pi_true = 3.14159265358979323846d0

  !======= Internals ============

  ! set the integrand function
  f(a) = 4.d0 / (1.d0 + a*a)

  ! intialize MPI
  call MPI_Init(ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Init =',ierr
     stop
  endif
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_size =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_rank =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! root outputs parallelism information and inputs the total number of intervals
  if (myid == 0) then
     print '(A,i5,A)', 'Running with ', numprocs,' MPI tasks'
     print *, 'Enter the total number of intervals (0 quits):'
     call flush()
     read(*,*) n
  endif

  ! root sends n out to other processors
  call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! stop for illegal n
  if (n < 1) then
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! start timer
  stime = MPI_Wtime()

  ! set subinterval width
  h = 1.d0/n

  ! determine this processor's interval
  is = floor(1.d0*n/numprocs)*myid + 1
  ie = floor(1.d0*n/numprocs)*(myid+1)
  if (myid == numprocs-1) then
     ie = n
  endif

  ! perform integration over n intervals
  mypi = 0.d0
  do i=is,ie,1
     x = h*(i - 0.5d0)
     mypi = mypi + h*f(x)
  enddo

  ! root node collects result
  call MPI_Reduce(mypi, pi, 1, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Reduce =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! stop timer
  ftime = MPI_Wtime();

  ! output computed value and error
  if (myid == 0) then
     print *, ' computed pi =',pi
     print *, '     true pi =',pi_true
     print *, '       error =',pi_true - pi
     print *, '     runtime =',ftime-stime
  endif

  ! finalize MPI
  call MPI_Finalize(ierr)

end program ComputePi
!=================================================================
