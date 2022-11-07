! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================

module exchange
  !===============================================================
  ! Description:
  !    Module containing data buffers for MPI exchanges with
  !    neighboring processes.  Also contains routines to:
  !      (a) initialize the module
  !      (b) perform exchanges
  !      (c) free module data.
  !===============================================================
  use mpi
  implicit none

  save       ! keep data available between use

  integer, private :: nxloc=0, nyloc=0
  integer, private :: nbN=-1, nbS=-1, nbE=-1, nbW=-1, comm=-1
  real*8, private, dimension(:), allocatable :: send_EW, send_NS

contains

  subroutine init_exchange(com, nxl, nyl, pcoords)
    !===============================================================
    ! Initializes module data
    !===============================================================
    ! inclusions
    implicit none

    ! declarations
    integer, intent(in) :: com, nxl, nyl, pcoords(2)
    integer :: nbcoords(2), ierr

    ! internals

    ! set communicator
    comm = com

    ! set subdomain extents
    nxloc = nxl
    nyloc = nyl

    ! allocate send buffers
    allocate(send_EW(nyloc), send_NS(nxloc))

    ! determine process neighbors
    nbcoords = (/ pcoords(1)-1, pcoords(2) /)
    call MPI_Cart_rank(comm, nbcoords, nbW, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'init_exchange error in MPI_Cart_rank =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    nbcoords = (/ pcoords(1)+1, pcoords(2) /)
    call MPI_Cart_rank(comm, nbcoords, nbE, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'init_exchange error in MPI_Cart_rank =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    nbcoords = (/ pcoords(1), pcoords(2)-1 /)
    call MPI_Cart_rank(comm, nbcoords, nbS, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'init_exchange error in MPI_Cart_rank =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    nbcoords = (/ pcoords(1), pcoords(2)+1 /)
    call MPI_Cart_rank(comm, nbcoords, nbN, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'init_exchange error in MPI_Cart_rank =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

  end subroutine init_exchange
  !===============================================================



  subroutine do_exchange1(v1, v1W, v1S, ierr)
    !===============================================================
    ! Performs neighbor exchange to update v1 buffers
    !===============================================================
    ! inclusions
    implicit none

    ! declarations
    integer, intent(out) :: ierr
    real*8, dimension(nxloc,nyloc), intent(in) :: v1
    real*8, intent(out) :: v1W(nyloc), v1S(nxloc)
    integer :: reqS, reqW, sreqN, sreqE
    integer :: status(MPI_STATUS_SIZE)

    ! internals

    ! check that module has been initialized
    if (nxloc < 1) then
       write(0,*) 'do_exchange1 error: exchange module not initialized!'
       call MPI_Abort(comm, 1, ierr)
    endif

    ! open receive channels for asynchronous communication
    call MPI_Irecv(v1W, nyloc, MPI_REAL8, nbW, 1, comm, reqW, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Irecv =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Irecv(v1S, nxloc, MPI_REAL8, nbS, 2, comm, reqS, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Irecv =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

    ! pack send buffers
    send_EW = v1(nxloc,1:nyloc)
    send_NS = v1(1:nxloc,nyloc)

    ! send data to neighbors
    call MPI_Isend(send_EW, nyloc, MPI_REAL8, nbE, 1, comm, sreqE, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Isend =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Isend(send_NS, nxloc, MPI_REAL8, nbN, 2, comm, sreqN, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Isend =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

    ! wait for communications to finish
    call MPI_Wait(sreqE, status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(sreqN, status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(reqW,  status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(reqS,  status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange1 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

  end subroutine do_exchange1
  !=================================================================



  subroutine do_exchange2(v2, v3, v2E, v3N, ierr)
    !===============================================================
    ! Performs neighbor exchange to update v2 and v3 buffers
    !===============================================================
    ! inclusions
    implicit none

    ! declarations
    integer, intent(out) :: ierr
    real*8, dimension(nxloc,nyloc), intent(in) :: v2, v3
    real*8, intent(out) :: v2E(nyloc), v3N(nxloc)
    integer :: reqN, reqS, reqE, reqW
    integer :: sreqN, sreqS, sreqE, sreqW
    integer :: status(MPI_STATUS_SIZE)

    ! internals

    ! check that module has been initialized
    if (nxloc < 1) then
       write(*,*) 'do_exchange2 error: exchange module not initialized!'
       return
    endif

    ! open receive channels for asynchronous communication
    call MPI_Irecv(v2E, nyloc, MPI_REAL8, nbE, 1, comm, reqE, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Irecv =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Irecv(v3N, nxloc, MPI_REAL8, nbN, 2, comm, reqN, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Irecv =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

    ! pack send buffers
    send_EW = v2(1,1:nyloc)
    send_NS = v3(1:nxloc,1)

    ! send data to neighbors
    call MPI_Isend(send_EW, nyloc, MPI_REAL8, nbW, 1, comm, sreqW, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Isend =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Isend(send_NS, nxloc, MPI_REAL8, nbS, 2, comm, sreqS, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Isend =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

    ! wait for communications to finish
    call MPI_Wait(sreqW, status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(sreqS, status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(reqE,  status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif
    call MPI_Wait(reqN,  status, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*) 'do_exchange2 error in MPI_Wait =',ierr
       call MPI_Abort(comm, 1, ierr)
    endif

  end subroutine do_exchange2
  !=================================================================



  subroutine free_exchange()
    !===============================================================
    ! Frees module data
    !===============================================================
    implicit none
    deallocate(send_EW, send_NS)

  end subroutine free_exchange
  !===============================================================

end module exchange
!===============================================================
