! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 9 January 2009
!=================================================================


program ChemicalEquilibrium_MPI
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the equilibrium chemical densities at a number of 
  !    spatial locations, given a (random) background temperature 
  !    field.  The chemical rate equations and solution strategy
  !    are in the subroutine chem_solver, which is called at every
  !    spatial location.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi
!  include "mpif.h"  ! include this if mpi module unavailable
  implicit none

  !======= Interfaces ===========
  interface
     subroutine chem_solver(T,u,v,w,lam,eps,maxit,its,res)
       integer, intent(in) :: maxit
       integer, intent(out) :: its
       double precision, intent(in) :: T, lam, eps
       double precision, intent(out) :: res
       double precision, intent(inout) :: u, v, w
     end subroutine chem_solver
  end interface

  !======= Declarations =========
  integer :: n, i, is, ie, js, je, maxit, its, numprocs, myid, ierr
  integer, allocatable :: counts(:), displs(:)
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: lam, eps, res, stime, ftime

  !======= Internals ============
  
  ! set solver input parameters
  maxit = 1000000
  lam = 1.d-2
  eps = 1.d-10
  
  ! intialize MPI
  call mpi_init(ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_init =',ierr
     stop
  endif
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_comm_size =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_comm_rank =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif

  ! input the total number of intervals
  if (myid == 0) then
     write(*,*) 'Enter the total number of intervals (0 quits):'
     read(*,*) n
  endif

  ! root node outputs parallelism information to screen
  if (myid == 0) then
     write(*,'(A,i5,A)') ' starting MPI with ', numprocs,' processes'
  endif

  ! root node sends n out to other processors
  call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_bcast =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif
  
  ! stop for illegal n
  if (n < 1) then
     call mpi_finalize(ierr)
     stop
  endif

  ! determine this processor's interval
  is = floor(1.d0*n/numprocs)*myid + 1
  ie = floor(1.d0*n/numprocs)*(myid+1)
  if (myid == numprocs-1) then
     ie = n
  endif

  ! root node allocates temperature field and solution arrays
  if (myid == 0) then
     allocate(T(n),u(n),v(n),w(n))

     ! set random temperature field, initial guesses at chemical densities
     call random_number(T)
     u = 0.35d0
     v = 0.1d0
     w = 0.5d0
  else  ! other nodes allocate local temperature field and solution arrays
     allocate(T(is:ie),u(is:ie),v(is:ie),w(is:ie))

     ! set initial guesses at chemical densities
     u = 0.35d0
     v = 0.1d0
     w = 0.5d0
  endif

  ! allocate gatherv/scatterv temporary arrays
  allocate(counts(numprocs),displs(numprocs))

  ! start timer
  stime = MPI_Wtime()

  ! root node sends out portions of temperature field to all procs
  do i=0,numprocs-1
     js = floor(1.d0*n/numprocs)*i + 1
     je = floor(1.d0*n/numprocs)*(i+1)
     if (i == numprocs-1) then
        je = n
     endif
     counts(i+1) = je-js+1
     displs(i+1) = js-1
  enddo
  call mpi_scatterv(T, counts, displs, MPI_DOUBLE_PRECISION, T, ie-is+1, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

  ! call solver over local intervals
  do i=is,ie,1
     call chem_solver(T(i),u(i),v(i),w(i),lam,eps,maxit,its,res)
     if (res < eps) then
!!$        write(*,'(2(A,i8))') '    i =',i,',  its =',its
     else
        write(*,'(2(A,i8),5(A,es9.2))') '    error: i =',i,',  its =',its, &
             ',  res =',res,',  T =',T(i),',  u =',u(i),',  v =',v(i),',  w =',w(i)
        call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     endif
  enddo

  ! root node collects results
  call mpi_gatherv(u, ie-is+1, MPI_DOUBLE_PRECISION, u, counts, displs, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_gatherv =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif
  call mpi_gatherv(v, ie-is+1, MPI_DOUBLE_PRECISION, v, counts, displs, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_gatherv =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif
  call mpi_gatherv(w, ie-is+1, MPI_DOUBLE_PRECISION, w, counts, displs, &
                   MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_gatherv =',ierr
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif

  ! stop timer
  ftime = MPI_Wtime()

  ! output solution time
  write(*,'(A,i3,A,i5,A)') ' proc ',myid,' computed ',ie-is+1,' iterations'
  if (myid == 0) then
     write(*,*) '     runtime =',ftime-stime
  endif


  ! free temperature/solution arrays and mpi temporary arrays
  deallocate(T,u,v,w,counts,displs)

  ! finalize MPI
  call mpi_finalize(ierr)  

end program ChemicalEquilibrium_MPI
!=================================================================
