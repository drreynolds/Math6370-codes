! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 10 March 2009
!=================================================================


program simple
  !-----------------------------------------------------------------
  ! Description: 
  !    This is a simple program using the basic 6 MPI functions.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi
!  include "mpif.h"  ! include this if mpi module unavailable

  !======= Declarations =========
  implicit none
  integer :: numprocs, myid, ierr, number, p
  integer :: status(MPI_STATUS_SIZE), tag, sender

  !======= Internals ============
  
  ! intialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)

  ! everyone send their (ID+1) to the root processor (except for root)
  !   - set the tag to mirror the sending processor id
  if (myid /= 0) then
     number = myid + 1
     call MPI_Send(number, 1, MPI_INTEGER, 0, myid, MPI_COMM_WORLD, ierr)
  endif

  ! the root node receives these (in order) and outputs each to screen
  if (myid == 0) then
     
     ! loop over all other processors
     do p=1,numprocs-1

        ! receive the number from this processor
        call MPI_Recv(number, 1, MPI_INTEGER, p, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, status, ierr)

        ! get some information about this message
        tag = status(MPI_TAG)
        sender = status(MPI_SOURCE)

        ! output the information and the data to screen
        print '(4(A,i4))', 'received value ',number,' from processor ', &
             p, ', tag =',tag,', sender =',sender

     end do
     
  endif
  
  ! finalize MPI
  call MPI_Finalize(ierr)  

end program simple
!=================================================================
