! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================


subroutine linresid2D(u, f, res, norm2, locN, locM, dx, dy, comm, ierr)
  !-----------------------------------------------------------------
  ! Description:
  !    calculates the 2D linear residual
  !             res = L*u - f
  !    and its L2-norm.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none

  integer, intent(in) :: locN, locM, comm
  integer, intent(out) :: ierr
  double precision, intent(in) :: dx, dy
  double precision, intent(in),  dimension(locN,locM) :: u, f
  double precision, intent(out), dimension(locN,locM) :: res
  double precision, intent(out) :: norm2
  double precision, dimension(locN) :: Nrecv, Srecv, Nsend, Ssend
  double precision, dimension(locM) :: Erecv, Wrecv, Esend, Wsend
  double precision :: tmp
  integer :: my_id, i, j, pdims(2), pcoords(2)
  logical :: periods(2)
  integer :: nbE, nbW, nbN, nbS, nbcoords(2)
  integer :: reqE, reqW, reqN, reqS, sreqE, sreqW, sreqN, sreqS
  integer :: status(MPI_STATUS_SIZE)

  !======= Internals ============

  ! Get MPI parallelism information from comm
  call MPI_Cart_get(comm, 2, pdims, periods, pcoords, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) ' error in MPI_Cart_get = ', ierr
     call MPI_Abort(comm, 1, ierr)
  endif
  call MPI_Comm_rank(comm, my_id, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) ' error in MPI_Comm_rank = ', ierr
     call MPI_Abort(comm, 1, ierr)
  endif

  ! initialize send/recv buffers
  Esend = u(locN,1:locM)
  Wsend = u(1,1:locM)
  Nsend = u(1:locN,locM)
  Ssend = u(1:locN,1)
  Nrecv = 0.d0
  Srecv = 0.d0
  Erecv = 0.d0
  Wrecv = 0.d0

  ! determine process 'neighbors'
  if (pcoords(1) /= 0) then
     nbcoords(1) = pcoords(1)-1
     nbcoords(2) = pcoords(2)
     call MPI_Cart_rank(comm, nbcoords, nbW, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Cart_rank = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  end if
  if (pcoords(1) /= pdims(1)-1) then
     nbcoords(1) = pcoords(1)+1
     nbcoords(2) = pcoords(2)
     call MPI_Cart_rank(comm, nbcoords, nbE, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Cart_rank = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  end if
  if (pcoords(2) /= 0) then
     nbcoords(1) = pcoords(1)
     nbcoords(2) = pcoords(2)-1
     call MPI_Cart_rank(comm, nbcoords, nbS, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Cart_rank = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  end if
  if (pcoords(2) /= pdims(2)-1) then
     nbcoords(1) = pcoords(1)
     nbcoords(2) = pcoords(2)+1
     call MPI_Cart_rank(comm, nbcoords, nbN, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Cart_rank = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  end if

  ! phase 1: open receive channels for neighbor values
  if (pcoords(1) /= 0) then
     call MPI_Irecv(Wrecv, locM, MPI_REAL8, nbW, 1, comm, reqW, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Irecv = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(1) /= pdims(1)-1) then
     call MPI_Irecv(Erecv, locM, MPI_REAL8, nbE, 2, comm, reqE, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Irecv = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= 0) then
     call MPI_Irecv(Srecv, locN, MPI_REAL8, nbS, 3, comm, reqS, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Irecv = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= pdims(2)-1) then
     call MPI_Irecv(Nrecv, locN, MPI_REAL8, nbN, 4, comm, reqN, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Irecv = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif

  ! phase 2: send boundary values to neighbors
  if (pcoords(1) /= pdims(1)-1) then
     call MPI_Isend(Esend, locM, MPI_REAL8, nbE, 1, comm, sreqE, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Isend = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(1) /= 0) then
     call MPI_Isend(Wsend, locM, MPI_REAL8, nbW, 2, comm, sreqW, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Isend = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= pdims(2)-1) then
     call MPI_Isend(Nsend, locN, MPI_REAL8, nbN, 3, comm, sreqN, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Isend = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= 0) then
     call MPI_Isend(Ssend, locN, MPI_REAL8, nbS, 4, comm, sreqS, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Isend = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif

  ! phase 3: wait until receives finish
  if (pcoords(1) /= 0) then
     call MPI_Wait(reqW, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(1) /= pdims(1)-1) then
     call MPI_Wait(reqE, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= 0) then
     call MPI_Wait(reqS, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= pdims(2)-1) then
     call MPI_Wait(reqN, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif

  ! compute linear residual in interior of subdomain
  norm2 = 0.d0
  do j=2,locM-1
     do i=2,locN-1
        res(i,j) = (u(i-1,j) - 2.d0*u(i,j) + u(i+1,j))/dx/dx &
             + (u(i,j-1) - 2.d0*u(i,j) + u(i,j+1))/dy/dy - f(i,j)
        norm2 = norm2 + dx*dy*res(i,j)**2
     end do
  end do

  ! compute linear residual at South edge of local subdomain
  do i=2,locN-1
     res(i,1) = (u(i-1,1) - 2.d0*u(i,1) + u(i+1,1))/dx/dx &
          + (Srecv(i) - 2.d0*u(i,1) + u(i,2))/dy/dy - f(i,1)
     norm2 = norm2 + dx*dy*res(i,1)**2
  end do

  ! compute linear residual at West edge of local subdomain
  do j=2,locM-1
     res(1,j) = (Wrecv(j) - 2.d0*u(1,j) + u(2,j))/dx/dx           &
          + (u(1,j-1) - 2.d0*u(1,j) + u(1,j+1))/dy/dy - f(1,j)
     norm2 = norm2 + dx*dy*res(1,j)**2
  end do

  ! compute linear residual at East edge of local subdomain
  do j=2,locM-1
     res(locN,j) = (u(locN-1,j) - 2.d0*u(locN,j) + Erecv(j))/dx/dx &
          + (u(locN,j-1) - 2.d0*u(locN,j) + u(locN,j+1))/dy/dy - f(locN,j)
     norm2 = norm2 + dx*dy*res(locN,j)**2
  end do

  ! compute linear residual at North edge of local subdomain
  do i=2,locN-1
     res(i,locM) = (u(i-1,locM) - 2.d0*u(i,locM) + u(i+1,locM))/dx/dx &
          + (u(i,locM-1) - 2.d0*u(i,locM) + Nrecv(i))/dy/dy - f(i,locM)
     norm2 = norm2 + dx*dy*res(i,locM)**2
  end do

  ! compute linear residual at corners of local subdomain
  res(1,1) = (Wrecv(1) - 2.d0*u(1,1) + u(2,1))/dx/dx &
           + (Srecv(1) - 2.d0*u(1,1) + u(1,2))/dy/dy - f(1,1)
  norm2 = norm2 + dx*dy*res(1,1)**2

  res(locN,1) = (u(locN-1,1) - 2.d0*u(locN,1) + Erecv(1))/dx/dx &
       + (Srecv(locN) - 2.d0*u(locN,1) + u(locN,2))/dy/dy - f(locN,1)
  norm2 = norm2 + dx*dy*res(locN,1)**2

  res(1,locM) = (Wrecv(locM) - 2.d0*u(1,locM) + u(2,locM))/dx/dx &
       + (u(1,locM-1) - 2.d0*u(1,locM) + Nrecv(1))/dy/dy - f(1,locM)
  norm2 = norm2 + dx*dy*res(1,locM)**2

  res(locN,locM) = (u(locN-1,locM) - 2.d0*u(locN,locM) + Erecv(locM))/dx/dx &
       + (u(locN,locM-1) - 2.d0*u(locN,locM) + Nrecv(locN))/dy/dy - f(locN,locM)
  norm2 = norm2 + dx*dy*res(locN,locM)**2


  ! combine local sums into global L2-norm
  call MPI_Allreduce(norm2, tmp, 1, MPI_REAL8, MPI_SUM, comm, ierr)
  if (ierr /= 0) then
     print *,'linresid error: MPI_Allreduce failed with error = ',ierr
     call MPI_Abort(comm, 1, ierr)
  end if
  norm2 = sqrt(tmp)


  ! last phase: wait until sends finish
  if (pcoords(1) /= 0) then
     call MPI_Wait(sreqW, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(1) /= pdims(1)-1) then
     call MPI_Wait(sreqE, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= 0) then
     call MPI_Wait(sreqS, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif
  if (pcoords(2) /= pdims(2)-1) then
     call MPI_Wait(sreqN, status, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) ' error in MPI_Wait = ', ierr
        call MPI_Abort(comm, 1, ierr)
     endif
  endif

  return

end subroutine linresid2D
!=================================================================
