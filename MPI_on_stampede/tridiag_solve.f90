! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 23 March 2009
!=================================================================


subroutine tridiag_solve(a, b, c, u, r, local_N, comm, ierr)
  !-----------------------------------------------------------------
  ! Description: 
  !   Solves the linear system 
  !           a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k),
  !   using a parallelized direct tridiagonal solver.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none

  integer, intent(in) :: local_N, comm
  integer, intent(out) :: ierr
  double precision, dimension(local_N) :: a, b, c, u, r
  double precision :: buff(3), b_l, r_l, c_l, u_r
  integer :: status(MPI_STATUS_SIZE)
  integer, parameter :: msg_xch_fwd = 0
  integer, parameter :: msg_xch_bck = 1
  integer :: k, my_id, nprocs

  !======= Internals ============
  
  ! get parallelism information from MPI communicator
  call MPI_Comm_size(comm, nprocs, ierr);
  call MPI_Comm_rank(comm, my_id, ierr);

  ! FORWARD SWEEP
  !   Wait for left b,r,c values from left neighbor processor
  if (my_id /= 0) then
     call MPI_Recv(buff, 3, MPI_DOUBLE_PRECISION, my_id-1, &
          msg_xch_fwd, comm, status, ierr)
     if (ierr /= 0) then
        print *,' tridiag_solve error: failed MPI_Recv with error = ',ierr
        return
     end if
     b_l = buff(1)
     r_l = buff(2)
     c_l = buff(3)
  else 
     b_l = 1.d0
     r_l = 0.d0
     c_l = 0.d0
  end if

  !   Perform local forward sweep, updating the diagonal and rhs
  b(1) = b(1) - c_l*a(1)/b_l;
  r(1) = r(1) - r_l*a(1)/b_l;
  do k=2,local_N
     b(k) = b(k) - c(k-1)*a(k)/b(k-1)
     r(k) = r(k) - r(k-1)*a(k)/b(k-1)
  end do

  !   Send right-most b,r,c values to right neighbor processor
  if (my_id /= nprocs-1) then
    buff(1) = b(local_N)
    buff(2) = r(local_N)
    buff(3) = c(local_N)
    call MPI_Send(buff, 3, MPI_DOUBLE_PRECISION, my_id+1, &
         msg_xch_fwd, comm, ierr)
    if (ierr /= 0) then
      print *, ' tridiag_solve error: failed MPI_Send with error = ',ierr
      return
    end if
  end if



  ! BACKWARD SWEEP
  !   Wait for right boundary value from right neighbor processor
  if (my_id /= nprocs-1) then
     call MPI_Recv(u_r, 1, MPI_DOUBLE_PRECISION, my_id+1, &
          msg_xch_bck, comm, status, ierr)
     if (ierr /= 0) then
        print *, ' tridiag_solve error: failed MPI_Recv with error = ',ierr
        return
     end if
  else 
     u_r = 0.d0
  end if

  !   Perform local backward sweep, updating the solution
  u(local_N) = (r(local_N) - c(local_N)*u_r)/b(local_N)
  do k=local_N-1,1,-1
     u(k) = (r(k)-c(k)*u(k+1))/b(k)
  end do
  
  !   Send left-most value to left neighbor processor
  if (my_id /= 0) then
     u_r = u(1)
     call MPI_Send(u_r, 1, MPI_DOUBLE_PRECISION, my_id-1, msg_xch_bck, &
          comm, ierr)
     if (ierr /= 0) then
        print *, ' tridiag_solve error: failed MPI_Send with error = ',ierr
        return
     end if
  end if

  ! Wait until all procs have caught up
  call MPI_Barrier(comm, ierr)
  if (ierr /= 0) then
     print *, ' tridiag_solve error: failed MPI_Barrier with error = ',ierr
     return
  end if

  ! finished
  return

end subroutine tridiag_solve
!=================================================================
