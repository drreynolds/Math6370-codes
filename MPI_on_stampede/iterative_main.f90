! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 23 March 2009
!=================================================================


program Iterative_main
  !-----------------------------------------------------------------
  ! Description: 
  !   We set up and solve the linear system 
  !           (I + gamma*L)u = r,
  !   where L is a standard 1D Laplace operator, r is a given 
  !   right-hand side, and u is the solution, using a parallelized 
  !   Jacobi iterative solver.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi
  implicit none

  !======= Interfaces ===========
  interface
     subroutine linresid(a, b, c, u, r, res, norm2, locN, N, comm, ierr)
       integer, intent(in) :: locN, N, comm
       integer, intent(out) :: ierr
       double precision, intent(in),  dimension(locN) :: a, b, c, u, r
       double precision, intent(out), dimension(locN) :: res
       double precision, intent(out) :: norm2
     end subroutine linresid
     subroutine jacobi_solve(a, b, c, u, r, res, local_N, global_N, &
          delta, maxiter, iters, comm, ierr)
       integer, intent(in) :: local_N, global_N, comm, maxiter
       integer, intent(out) :: ierr
       double precision, intent(in), dimension(local_N) :: a, b, c, r
       double precision, dimension(local_N) :: u, res
       double precision, intent(in) :: delta
       integer, intent(out) :: iters
     end subroutine jacobi_solve
  end interface

  !======= Declarations =========
  double precision, allocatable :: u(:), r(:), a(:), b(:), c(:), res(:)
  double precision :: gamma, delta
  integer :: local_N, global_N, k, my_id, nprocs, ierr, comm
  integer :: maxiter, iters
  double precision :: stime, ftime, err2norm
  namelist /inputs/ gamma, delta, global_N

  !======= Internals ============
  
  ! intialize MPI
  call mpi_init(ierr)
  if (ierr /= 0) then
     write(*,*) ' error in mpi_init =',ierr
     stop
  endif
  comm = MPI_COMM_WORLD
  call mpi_comm_size(comm, nprocs, ierr)
  call mpi_comm_rank(comm, my_id, ierr)

  ! root node gets problem information from input file
  if (my_id == 0) then
     open(100,file='input_iterative.txt',form='formatted')
     read(100,inputs)
     close(100)
  end if

  ! root node broadcasts information to other procs
  call mpi_bcast(gamma, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  if (ierr /= 0) then
     print *, ' error in mpi_bcast =',ierr
     call MPI_Abort(comm, -1, ierr);
     stop
  endif
  call mpi_bcast(delta, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr)
  if (ierr /= 0) then
     print *, ' error in mpi_bcast =',ierr
     call MPI_Abort(comm, -1, ierr);
     stop
  endif
  call mpi_bcast(global_N, 1, MPI_INTEGER, 0, comm, ierr)
  if (ierr /= 0) then
     print *, ' error in mpi_bcast =',ierr
     call MPI_Abort(comm, -1, ierr);
     stop
  endif

  ! set local mesh sizes (last process takes up remainder)
  local_N = global_N/nprocs
  if (my_id == nprocs-1) then
     if (local_N*nprocs /= global_N) then
        local_N = global_N - local_N*(nprocs-1)
     end if
  end if

  ! root node outputs some information to screen
  if (my_id == 0) then
       print *,'iterative test with ',nprocs,' processors'
       print *,'    gamma = ',gamma
       print *,'    linear solver tolerance delta = ',delta
       print *,'    global problem size N = ',global_N
       print *,'    local problem sizes n = ',local_N
  endif

  ! Allocate vector memory
  allocate(u(local_N), r(local_N), a(local_N), b(local_N), &
           c(local_N), res(local_N), stat=ierr)
  if (ierr /= 0) then
     print *,'iterative test error: proc ',my_id,&
          ', failed to allocate array data'
     call MPI_Abort(comm, -1, ierr)
     stop
  end if

  ! Set up matrix arrays, right-hand side, and initial solution guess
  do k=1,local_N
     u(k) =  0.d0
     r(k) =  1.d0
     a(k) = -gamma
     b(k) =  1.d0+gamma*2.d0
     c(k) = -gamma
  end do
  
  ! Adjust a, c arrays if we are at either end of domain
  if (my_id == 0) then
     a(1) = 0.d0
  end if
  if (my_id == nprocs-1) then
     c(local_N) = 0.d0
  end if


  ! check linear residual
  call linresid(a, b, c, u, r, res, err2norm, local_N, global_N, comm, ierr)
  if (ierr /= 0) then
     print *,'iterative test error: linresid failed with error = ',ierr
     call MPI_Abort(comm, -1, ierr)
     stop
  end if
  if (my_id == 0) then
     print *,' initial residual: ||T*u-r||_2 = ',err2norm
  end if

  ! Wait until all procs have caught up
  call MPI_Barrier(comm, ierr)
  if (ierr /= 0) then
     print *,'direct test error: failed MPI_Barrier with error = ',ierr
     call MPI_Abort(comm, -1, ierr)
     stop
  end if

  ! Solve system, get timing information
  maxiter = 10000
  stime = MPI_Wtime()
  call jacobi_solve(a, b, c, u, r, res, local_N, global_N, delta, &
       maxiter, iters, comm, ierr)
  if (ierr /= 0) then
     print *,'iterative test error: jacobi_solve failed with error = ',ierr
     call MPI_Abort(comm, -1, ierr)
     stop
  end if
  ftime = MPI_Wtime()
  if (my_id == 0) then
     print *,' converged in ',iters,' iterations at delta = ',delta
     print *,' solution time: ',ftime-stime,' seconds'
  end if
  

  ! check linear residual
  call linresid(a, b, c, u, r, res, err2norm, local_N, global_N, comm, ierr)
  if (ierr /= 0) then
     print *,'iterative test error: linresid failed with error = ',ierr
     call MPI_Abort(comm, -1, ierr)
     stop
  end if
  if (my_id == 0) then
     print *,' final residual: ||T*u-r||_2 = ',err2norm
  end if

  
  ! Free matrix/solver memory
  deallocate(u,r,a,b,c,res)

  ! finalize MPI
  call MPI_Finalize(ierr)

end program Iterative_main
!=================================================================
