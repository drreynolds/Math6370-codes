! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 23 March 2009
!=================================================================


subroutine linresid(a, b, c, u, r, res, norm2, locN, N, comm, ierr)
  !-----------------------------------------------------------------
  ! Description: 
  !    calculates the linear residual and its averaged 2-norm (WRMS)
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none

  integer, intent(in) :: locN, N, comm
  integer, intent(out) :: ierr
  double precision, intent(in),  dimension(locN) :: a, b, c, u, r
  double precision, intent(out), dimension(locN) :: res
  double precision, intent(out) :: norm2
  double precision :: u_l, u_r, s_l, s_r
  double precision :: tmp
  integer :: my_id, nprocs, k
  integer :: status(MPI_STATUS_SIZE)
  integer :: id_recv_l, id_recv_r, id_send_l, id_send_r

  !======= Internals ============
  
  ! Get MPI parallelism information from comm
  call MPI_Comm_size(comm, nprocs, ierr)
  call MPI_Comm_rank(comm, my_id, ierr)

  ! open asynchronous receive channels for neighbor boundary values
  u_l = 0.d0
  u_r = 0.d0
  if (my_id /= 0) then
     call MPI_Irecv(u_l, 1, MPI_DOUBLE_PRECISION, my_id-1, 100, &
          comm, id_recv_l, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Irecv failed with error = ',ierr
        return
     end if
  end if
  if (my_id /= nprocs-1) then
     call MPI_Irecv(u_r, 1, MPI_DOUBLE_PRECISION, my_id+1, 101, &
          comm, id_recv_r, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Irecv failed with error = ',ierr
        return
     end if
  end if

  ! send boundary values to neighbor processes
  s_l = u(1)
  s_r = u(locN)
  if (my_id /= nprocs-1) then
     call MPI_Isend(s_r, 1, MPI_DOUBLE_PRECISION, my_id+1, 100, &
          comm, id_send_r, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Isend failed with error = ',ierr
        return
     end if
  end if
  if (my_id /= 0) then
     call MPI_Isend(s_l, 1, MPI_DOUBLE_PRECISION, my_id-1, 101, &
          comm, id_send_l, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Isend failed with error = ',ierr
        return
     end if
  end if

  ! wait for left boundary value to arrive before using
  if (my_id /= 0) then
     call MPI_Wait(id_recv_l, status, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Wait failed with error = ',ierr
        return
     end if
  end if

  ! compute linear residual at left of local subdomain
  res(1) = a(1)*u_l + b(1)*u(1) + c(1)*u(2) - r(1)
  norm2 = res(1)**2

  ! compute linear residual in interior of subdomain
  do k=2,locN-1
    res(k) = a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) - r(k)
    norm2 = norm2 + res(k)**2
  end do

  ! wait for right boundary value to arrive before using
  if (my_id /= nprocs-1) then
     call MPI_Wait(id_recv_r, status, ierr)
     if (ierr /= 0) then
        print *,'linresid error: MPI_Wait failed with error = ',ierr
        return
     end if
  end if

  ! compute linear residual at right end of local subdomain
  k = locN
  res(k) = a(k)*u(k-1) + b(k)*u(k) + c(k)*u_r - r(k)
  norm2 = norm2 + res(k)**2

  ! combine local 2-norms into global averaged 2-norm
  call MPI_Allreduce(norm2, tmp, 1, MPI_DOUBLE_PRECISION, &
       MPI_SUM, comm, ierr)
  if (ierr /= 0) then
     print *,'linresid error: MPI_Allreduce failed with error = ',ierr
     return
  end if
  norm2 = sqrt(tmp/N)

  return

end subroutine linresid
!=================================================================
