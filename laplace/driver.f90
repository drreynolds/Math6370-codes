! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 30 March 2009
!=================================================================


program Laplace2D
  !-----------------------------------------------------------------
  ! Description: 
  !   We set up the linear residual 
  !           L*u - f,
  !   where L is a standard 2D Laplace operator, and f and u are 
  !   given vectors.  These are decomposed using a 2D parallel 
  !   domain decomposition strategy.  We then call the routine 
  !   linresid2D to compute the linear residual above.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi
  implicit none

  !======= Interfaces ===========
  interface
     subroutine linresid2D(u,f,res,norm2,locN,locM,dx,dy,comm,ierr)
       integer, intent(in) :: locN, locM, comm
       integer, intent(out) :: ierr
       double precision, intent(in) :: dx, dy
       double precision, intent(in),  dimension(locN,locM) :: u, f
       double precision, intent(out), dimension(locN,locM) :: res
       double precision, intent(out) :: norm2
     end subroutine linresid2D
  end interface

  !======= Declarations =========
  double precision, allocatable :: u(:,:), f(:,:), res(:,:)
  integer :: N, M, locN, locM, i, j
  integer :: my_id, nprocs, px, py, ierr, comm
  integer :: buf(4), pdims(2), pcoords(2)
  logical :: periods(2)
  double precision :: dx, dy, xl, yl, x, y
  double precision :: stime, ftime, res2norm
  character*50 :: outfile
  namelist /inputs/ N, M, px, py

  !======= Internals ============

  ! intialize MPI
  call MPI_Init(ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Init =',ierr
     stop
  endif
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_size =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, my_id, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_rank =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! root gets problem information from input file
  if (my_id == 0) then
     open(100,file='input.txt',form='formatted')
     read(100,inputs)
     close(100)
     buf = (/ N, M, px, py /)
  end if

  ! root broadcasts information to other procs
  call MPI_Bcast(buf, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     print *, ' error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! extract values from buffer
  N  = buf(1)
  M  = buf(2)
  px = buf(3)
  py = buf(4)

  ! check that the processor layout and communicator size agree
  if (nprocs /= px*py) then
     if (my_id == 0) then
        write(0,*) ' error: incorrect processor layout'
        write(0,*) '     nprocs = ',nprocs
        write(0,*) '     px = ',px,',  py =',py
     endif
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! set up Cartesian communicator, comm
  pdims = (/ px, py /)
  periods = .false.
  call MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, .false., comm, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Cart_create =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! get this processor's new id and location from comm
  call MPI_Comm_rank(comm, my_id, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_rank =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Cart_coords(comm, my_id, 2, pcoords, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Cart_coords =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! set local mesh sizes (last process in each dir takes up remainder)
  dx = 1.d0/N
  dy = 1.d0/M
  locN = N/px
  locM = M/py
  xl = dx*locN*pcoords(1)
  yl = dy*locM*pcoords(2)
  if (pcoords(1) == px-1) then
     if (locN*px /= N) then
        locN = N - locN*(px-1)
     end if
  end if
  if (pcoords(2) == py-1) then
     if (locM*py /= M) then
        locM = M - locM*(py-1)
     end if
  end if

  ! root node outputs some information to screen
  if (my_id == 0) then
       print *,'Laplace2D driver with ',nprocs,' processes'
       print *,'    global problem size = ',N,'  by ',M
       print *,'    local problem sizes = ',locN,'  by ',locM
  endif

  ! Allocate vector memory
  allocate(u(locN,locM), f(locN,locM), res(locN,locM), stat=ierr)
  if (ierr /= 0) then
     write(0,*) 'Laplace2D error: proc ',my_id, &
          ', failed to allocate array data'
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  ! start timer
  stime = MPI_Wtime()

  ! Set up u and f
  do j=1,locM
     y = yl + (j-1)*dy
     do i=1,locN
        x = xl + (i-1)*dx
        u(i,j) = exp(-2.d1*((x-0.5d0)**2 + (y-0.5d0)**2))
        f(i,j) = 1.d0
     end do
  end do
  
  ! Adjust u at domain boundary
  if (pcoords(1) == 0) then
     u(1,:) = 0.d0
  end if
  if (pcoords(1) == px-1) then
     u(locN,:) = 0.d0
  end if
  if (pcoords(2) == 0) then
     u(:,1) = 0.d0
  end if
  if (pcoords(2) == py-1) then
     u(:,locM) = 0.d0
  end if

  ! Wait until all procs have caught up
  call MPI_Barrier(comm, ierr)
  if (ierr /= 0) then
     write(0,*) 'Laplace2D error: failed MPI_Barrier with error = ',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if

  ! Compute linear residual
  call linresid2D(u, f, res, res2norm, locN, locM, dx, dy, comm, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in linresid2D = ',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  if (my_id == 0) then
     print *,' residual: ||L*u-f||_2 = ',res2norm
  end if

  ! stop timer
  ftime = MPI_Wtime()
  if (my_id == 0) then
     print *,' runtime = ',ftime-stime
  end if

  
  ! output residual to file(s)
  !   root outputs problem/parallelism information
  if (my_id == 0) then
     open(101,file='resid.txt')
     write(101,*) N, M, px, py
     close(101)
  endif
  write(outfile,'(5Hresid,f4.3)') float(my_id)/1000.0
  open(102,file=outfile)
  write(102,*) locN
  write(102,*) locM
  write(102,*) pcoords(1)
  write(102,*) pcoords(2)
  do j=1,locM
     do i=1,locN
        write(102,*) res(i,j)
     end do
  end do
  close(102)
  
  ! Free allocated memory
  deallocate(u,f,res)

  ! finalize MPI
  call MPI_Finalize(ierr)

end program Laplace2D
!=================================================================
