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
  integer, parameter :: chunk=10
  integer :: n, i, ibase, maxit, its, numprocs, myid, ierr, mycount
  integer :: numsent, sender, ansentry, status(MPI_STATUS_SIZE)
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: buffer(3*chunk)
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

     ! do master/slave computation
     numsent = 0

     ! first send a part to each slave
     do i = 1,numprocs-1
        ! fill send buffer 
        if (numsent + chunk <= n) then
           buffer(1:chunk) = T(numsent+1:numsent+chunk)
        else
           buffer(1:n-numsent) = T(numsent+1:n)
           buffer(n-numsent+1:chunk) = 1.d0
        endif
        ! send with tag as entry
        call mpi_send(buffer, chunk, MPI_DOUBLE_PRECISION, i, i, &
                      MPI_COMM_WORLD, ierr)
        if (ierr /= 0) then
           write(*,*) ' error in mpi_send =',ierr
           call mpi_abort(MPI_COMM_WORLD, -1, ierr)
           stop
        endif
        numsent = numsent + chunk
     enddo

     ! if all work was done in first phase, receive work only
     if (numsent >= n) then

        ! receive work from each proc
        do i = 1,numprocs-1

           ! receive answers from any processor
           call mpi_recv(buffer, 3*chunk, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
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
           if (ansentry+chunk-1 <= n) then
              u(ansentry:ansentry+chunk-1) = buffer(1:chunk)
              v(ansentry:ansentry+chunk-1) = buffer(chunk+1:2*chunk)
              w(ansentry:ansentry+chunk-1) = buffer(2*chunk+1:3*chunk)
           else
              u(ansentry:n) = buffer(1:n-ansentry+1)
              v(ansentry:n) = buffer(chunk+1:n-ansentry+chunk+1)
              w(ansentry:n) = buffer(2*chunk+1:n-ansentry+2*chunk+1)
           endif

           ! tell senders that work is complete
           call mpi_send(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, &
                         0, MPI_COMM_WORLD, ierr)
           if (ierr /= 0) then
              write(*,*) ' error in mpi_send =',ierr
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              stop
           endif
           
        enddo

     ! receive work from each proc, and assign additional work if available
     else
        do i = 1,ceiling(1.d0*n/chunk)

           ! receive answers from any processor
           call mpi_recv(buffer, 3*chunk, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
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
           if (ansentry+chunk-1 <= n) then
              u(ansentry:ansentry+chunk-1) = buffer(1:chunk)
              v(ansentry:ansentry+chunk-1) = buffer(chunk+1:2*chunk)
              w(ansentry:ansentry+chunk-1) = buffer(2*chunk+1:3*chunk)
           else
              u(ansentry:n) = buffer(1:n-ansentry+1)
              v(ansentry:n) = buffer(chunk+1:n-ansentry+chunk+1)
              w(ansentry:n) = buffer(2*chunk+1:n-ansentry+2*chunk+1)
           endif
           if (numsent < n) then  ! send another set of rows
              if (numsent + chunk <= n) then
                 buffer(1:chunk) = T(numsent+1:numsent+chunk)
              else
                 buffer(1:n-numsent) = T(numsent+1:n)
                 buffer(n-numsent+1:chunk) = 1.d0
              endif
              call mpi_send(buffer, chunk, MPI_DOUBLE_PRECISION, sender, &
                            numsent+1, MPI_COMM_WORLD, ierr)
              if (ierr /= 0) then
                 write(*,*) ' error in mpi_send =',ierr
                 call mpi_abort(MPI_COMM_WORLD, -1, ierr)
                 stop
              endif
              numsent = numsent + chunk
           else   ! tell senders that work is complete
              call mpi_send(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION, sender, &
                            0, MPI_COMM_WORLD, ierr)
              if (ierr /= 0) then
                 write(*,*) ' error in mpi_send =',ierr
                 call mpi_abort(MPI_COMM_WORLD, -1, ierr)
                 stop
              endif
           endif
           
        enddo

     endif

     ! stop timer
     ftime = MPI_Wtime()


  else  ! slave role

     ! allocate temperature field and solution arrays
     allocate(T(chunk),u(chunk),v(chunk),w(chunk))

     ! loop over assignments from master
     more_work = .true.
     do while(more_work) 

        ! receive message from master
        call mpi_recv(buffer, chunk, MPI_DOUBLE_PRECISION, 0, MPI_ANY_TAG, &
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
           ibase = status(MPI_TAG)

           ! set input and initial guess
           T = buffer(1:chunk)
           u = 0.35d0
           v = 0.1d0
           w = 0.5d0
           
           ! call solver over entries
           do i=1,chunk
              call chem_solver(T(i),u(i),v(i),w(i),lam,eps,maxit,its,res)
              if (res < eps) then
!!$                 write(*,'(2(A,i8))') '    i =',ibase+i-1,',  its =',its
              else
                 write(*,'(2(A,i8),5(A,es9.2))') '    error: i =',ibase+i-1, &
                      ',  its =',its,',  res =',res,',  T =',T(i),',  u =',u(i), &
                      ',  v =',v(i),',  w =',w(i)
                 call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              endif
           enddo

           ! send result to master
           buffer(1:chunk) = u(1:chunk)
           buffer(chunk+1:2*chunk) = v(1:chunk)
           buffer(2*chunk+1:3*chunk) = w(1:chunk)
           call mpi_send(buffer, 3*chunk, MPI_DOUBLE_PRECISION, 0, ibase, &
                         MPI_COMM_WORLD, ierr)
           if (ierr /= 0) then
              write(*,*) ' error in mpi_send =',ierr
              call mpi_abort(MPI_COMM_WORLD, -1, ierr)
              stop
           endif
           
           ! update work counter
           mycount = mycount + chunk

        end if

     end do

  endif

  ! output solution time
  write(*,'(A,i3,A,i5,A)') ' proc ',myid,' computed ',mycount,' iterations'
  if (myid == 0) then
     write(*,*) '     runtime =',ftime-stime
  endif


  ! free temperature/solution arrays and mpi temporary arrays
  deallocate(T,u,v,w)

  ! finalize MPI
  call mpi_finalize(ierr)  

end program ChemicalEquilibrium_MPI
!=================================================================
