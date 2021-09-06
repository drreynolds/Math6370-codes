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
  integer :: n, i, maxit, its, numprocs, myid, ierr, mycount
  integer :: numsent, sender, ansentry, status(MPI_STATUS_SIZE)
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: Tval, uval, vval, wval
  double precision :: buffer(3)
  double precision :: lam, eps, res, stime, ftime
  logical :: more_work
  
  !======= Internals ============
  
  ! set solver input parameters
  maxit = 1000000
  lam = 1.d-2
  eps = 1.d-10

  ! initialize internal work counter
  mycount = 0
  
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

  ! check whether we have enough processors
  if (numprocs < 2) then
     write(*,*) ' master/slave example, requires > 1 proc!  Exiting'
     call mpi_abort(MPI_COMM_WORLD, -1, ierr)
     stop
  endif

  ! master role
  if (myid == 0) then

     ! input the total number of intervals
     write(*,*) 'Enter the total number of intervals (0 quits):'
     read(*,*) n
     write(*,'(A,i5,A)') ' starting MPI with ', numprocs,' processes'

     ! stop for illegal n
     if (n < 1) then
        call mpi_abort(MPI_COMM_WORLD, -1, ierr)
        stop
     endif

     ! allocate temperature field and solution arrays
     allocate(T(n),u(n),v(n),w(n))

     ! set random temperature field, initial guesses at chemical densities
     call random_number(T)

     ! start timer
     stime = MPI_Wtime()

     ! do master/slave computation (assumes n >= numprocs)
     numsent = 0

     ! first send a part to each slave
     do i = 1,min(n,numprocs-1)
        ! fill send buffer 
        buffer(1) = T(i)
        ! send with tag as entry
        call mpi_send(buffer, 1, MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) then
           write(*,*) ' error in mpi_send =',ierr
           call mpi_abort(MPI_COMM_WORLD, -1, ierr)
           stop
        endif
        numsent = numsent + 1
     enddo

     ! receive work from each proc, and assign additional work if available
     do i = 1,n

        ! receive answers from any processor
        call mpi_recv(buffer, 3, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) then
           write(*,*) ' error in mpi_recv =',ierr
           call mpi_abort(MPI_COMM_WORLD, -1, ierr)
           stop
        endif
        ! decode the sender and solution entry from status
        sender = status(MPI_SOURCE)
        ansentry = status(MPI_TAG)
        ! store results
        u(ansentry) = buffer(1)
        v(ansentry) = buffer(2)
        w(ansentry) = buffer(3)
        if (numsent < n) then  ! send another row
           buffer(1) = T(numsent+1)
           call mpi_send(buffer, 1, MPI_DOUBLE_PRECISION, sender, &
                         numsent+1, MPI_COMM_WORLD, ierr)
           if (ierr /= 0) then
              write(*,*) ' error in mpi_send =',ierr
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              stop
           endif
           numsent = numsent + 1
        else   ! tell senders that work is complete
           call mpi_send(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, 0, &
                         MPI_COMM_WORLD, ierr)
           if (ierr /= 0) then
              write(*,*) ' error in mpi_send =',ierr
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              stop
           endif
        endif
           
     enddo

     ! stop timer
     ftime = MPI_Wtime()

     ! free temperature/solution arrays
     deallocate(T,u,v,w)

  else  ! worker role

     ! loop over assignments from master
     more_work = .true.
     do while(more_work) 

        ! receive message from master
        call mpi_recv(buffer, 1, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, status, ierr)
        if (ierr /= 0) then
           write(*,*) ' error in mpi_recv =',ierr
           call mpi_abort(MPI_COMM_WORLD, -1, ierr)
           stop
        endif

        ! check if work is complete
        if (status(MPI_TAG) == 0) then
           more_work = .false.

        ! otherwise, do computation
        else

           ! set row index
           i = status(MPI_TAG)

           ! set input and initial guess
           Tval = buffer(1)
           uval = 0.35d0
           vval = 0.1d0
           wval = 0.5d0
           
           ! call solver
           call chem_solver(Tval,uval,vval,wval,lam,eps,maxit,its,res)
           if (res < eps) then
!!$              write(*,'(2(A,i8))') '    i =',i,',  its =',its
           else
              write(*,'(2(A,i8),5(A,es9.2))') '    error: i =',i,',  its =',its, &
                   ',  res =',res,',  T =',Tval,',  u =',uval,',  v =',vval,',  w =',wval
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
           endif
           
           ! send result to master
           buffer(1) = uval
           buffer(2) = vval
           buffer(3) = wval
           call mpi_send(buffer, 3, MPI_DOUBLE_PRECISION, 0, i, MPI_COMM_WORLD, ierr)
           if (ierr /= 0) then
              write(*,*) ' error in mpi_send =',ierr
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              stop
           endif
           
           ! update work counter
           mycount = mycount + 1

        endif

     end do

  endif

  ! output solution time
  write(*,'(A,i3,A,i5,A)') ' proc ',myid,' computed ',mycount,' iterations'
  if (myid == 0) then
     write(*,*) '     runtime =',ftime-stime
  endif


  ! finalize MPI
  call mpi_finalize(ierr)  

end program ChemicalEquilibrium_MPI
!=================================================================
