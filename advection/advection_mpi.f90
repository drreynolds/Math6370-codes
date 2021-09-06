! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!-----------------------------------------------------------------
! Description: 
!    Evolves the first-order 2D wave equations in time.
!=================================================================

program advection_mpi
  !===============================================================
  ! inclusions
  use mpi
  use exchange
  implicit none

  ! interfaces
  interface
     subroutine output_mpi(u, t, nxloc, nyloc, pid, pcoords, &
                           nx, ny, px, py, noutput)
       integer, intent(in) :: nxloc, nyloc, pid, pcoords(2)
       integer, intent(in) :: nx, ny, px, py, noutput
       real*8,  intent(in) :: t, u(nxloc,nyloc)
     end subroutine output_mpi
     subroutine initialize_mpi(u,v1,v2,v3,c,dx,dy,xl,yl,nxloc,nyloc,myid)
       integer, intent(in)  :: nxloc, nyloc, myid
       real*8,  intent(in)  :: c, dx, dy, xl, yl
       real*8, dimension(nxloc,nyloc), intent(out) :: u, v1, v2, v3
     end subroutine initialize_mpi
  end interface

  ! declarations
  integer :: nx, ny, nt, j, i, noutput, it, ierr
  integer :: myid, nprocs, ibuff(3), pdims(2), comm
  integer :: pcoords(2), nxloc, nyloc
  real*8  :: tstop, c, dtoutput, dx, dy, y_c, y_h, x_c, x_h, t, toutput, dt
  real*8  :: v2_E, v2_W, v3_N, v3_S, v1_W, v1_E, v1_S, v1_N
  real*8  :: stime, ftime, runtime, iotime, inittime, commtime
  real*8, dimension(:,:), allocatable :: u, v1, v2, v3
  real*8  :: rbuff(3), xl, yl 
  real*8, dimension(:), allocatable :: v2E, v3N, v1W, v1S
  logical :: periods(2)
  namelist /inputs/ nx, ny, nt, tstop, c, dtoutput

  ! internals
  call MPI_Init(ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Init =',ierr
     stop
  endif
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Comm_size =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Comm_rank =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! initialize timers
  runtime = 0.d0
  iotime = 0.d0
  inittime = 0.d0 
  commtime = 0.d0
  
  stime = MPI_Wtime()
  ! read problem parameters from input file (should be in this order):
  !    nx - number of nodes in x-direction
  !    ny - number of nodes in y-direction
  !    nt - number of time steps
  !    tstop - final time (will stop at nt or stop, whichever is 1st)
  !    c - wave speed
  !    dtoutput - time frequency for outputting solutions
  if (myid == 0) then
     open(100,file='input.txt')
     read(100,inputs)
     close(100)
     ibuff = (/ nx, ny, nt /)
     rbuff = (/ tstop, c, dtoutput /)
  end if
  call MPI_Bcast(ibuff, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Bcast(rbuff, 3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  nx = ibuff(1)
  ny = ibuff(2)
  nt = ibuff(3)
  tstop = rbuff(1)
  c = rbuff(2)
  dtoutput = rbuff(3)
  
  if (myid == 0) then
     print *, '  '
     print *, 'Running wave problem:'
     print *, '   nx =',nx,',  ny =',ny
     print *, '   nt =',nt,',  tstop =',tstop
     print *, '   c =',c
     print *, '   dtoutput =',dtoutput
     print *, '   nprocs =',nprocs
     print *, '  '
  end if
  ftime = MPI_Wtime()
  iotime = iotime + ftime-stime

  stime = MPI_Wtime()
  ! create cartesian communicator
  pdims = (/ 0, 0 /)
  periods = .true.
  call MPI_Dims_create(nprocs, 2, pdims, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Dims_create =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, &
       .false., comm, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Cart_create =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Cart_coords(comm, myid, 2, pcoords, ierr)
  if (ierr /= MPI_SUCCESS) then
     write(0,*) 'advection_mpi error in MPI_Cart_coords =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! set grid spacing and mesh points
  dx = 1.d0/nx
  dy = 1.d0/ny
  nxloc = nx/pdims(1)
  nyloc = ny/pdims(2)
  xl = dx*nxloc*pcoords(1)
  yl = dy*nyloc*pcoords(2)
  if (pcoords(1) == pdims(1)-1) then
     if (nxloc*pdims(1) /= nx) then
        nxloc = nx - nxloc*(pdims(1)-1)
     end if
  end if
  if (pcoords(2) == pdims(2)-1) then
     if (nyloc*pdims(2) /= ny) then
        nyloc = ny - nyloc*(pdims(2)-1)
     end if
  end if

  ! allocate arrays
  allocate(u(nxloc,nyloc), v1(nxloc,nyloc), &
           v2(nxloc,nyloc), v3(nxloc,nyloc))
  allocate(v2E(nyloc), v3N(nxloc), v1W(nyloc), v1S(nxloc))

  ! set initial conditions
  call initialize_mpi(u,v1,v2,v3,c,dx,dy,xl,yl,nxloc,nyloc,myid)
  ftime = MPI_Wtime()
  inittime = inittime + ftime-stime

  ! set initial time, output initial solution
  !   root outputs problem/parallelism information
  stime = MPI_Wtime()
  t = 0.d0
  toutput = 0.d0
  noutput = 0
  if (myid == 0) print '(A,i4,A,i8,A,es9.2)', &
       'writing output file ',noutput, ', step = ',0,',  t =',t
  call output_mpi(u, t, nxloc, nyloc, myid, pcoords, nx, ny, pdims(1), pdims(2), noutput)
  ftime = MPI_Wtime()
  iotime = iotime + ftime-stime

  ! initialize exchange module
  stime = MPI_Wtime()
  call init_exchange(comm, nxloc, nyloc, pcoords)
  ftime = MPI_Wtime()
  inittime = inittime + ftime-stime

  ! set time step
  dt = min(dx/c/5.d1, dy/c/5.d1)

  ! start time stepping
  do it=1,nt

     ! communicate to get v2 and v3 neighbor values
     stime = MPI_Wtime()
     call do_exchange2(v2, v3, v2E, v3N, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) 'advection_mpi error in do_exchange2 =',ierr
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     endif
     ftime = MPI_Wtime()
     commtime = commtime + ftime-stime

     stime = MPI_Wtime()
     ! first update v1 to get to half time step
     !    get relevant values for this location
     do j=1,nyloc
        do i=1,nxloc

           ! access relevant components of v2 and v3
           if (i==nxloc) then
              v2_E = v2E(j)
           else
              v2_E = v2(i+1,j)
           end if
           v2_W = v2(i,j)
           if (j==nyloc) then
              v3_N = v3N(i)
           else
              v3_N = v3(i,j+1)
           end if
           v3_S = v3(i,j)

           ! update v1
           v1(i,j) = v1(i,j) + c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S)

        end do
     end do
     ftime = MPI_Wtime()
     runtime = runtime + ftime-stime

     ! communicate to get v1 neighbor values
     stime = MPI_Wtime()
     call do_exchange1(v1, v1W, v1S, ierr)
     if (ierr /= MPI_SUCCESS) then
        write(0,*) 'advection_mpi error in do_exchange1 =',ierr
        call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
     endif
     ftime = MPI_Wtime()
     commtime = commtime + ftime-stime

     stime = MPI_Wtime()
     ! next update v2 & v3 to get to full time step
     !    get relevant values for this location
     do j=1,nyloc
        do i=1,nxloc

           ! access relevant components of v1
           if (i==1) then
              v1_W = v1W(j)
           else
              v1_W = v1(i-1,j)
           end if
           v1_E = v1(i,j)
           if (j==1) then
              v1_S = v1S(i)
           else
              v1_S = v1(i,j-1)
           end if
           v1_N = v1(i,j)

           ! update v2 and v3
           v2(i,j) = v2(i,j) + c*dt/dx*(v1_E - v1_W)
           v3(i,j) = v3(i,j) + c*dt/dy*(v1_N - v1_S)

        end do
     end do

     ! update solution
     u = u + dt*v1

     ! update time
     t = t + dt
     ftime = MPI_Wtime()
     runtime = runtime + ftime-stime

     ! stop simulation if we've reached tstop
     if (t >= tstop)  exit

     ! output solution periodically
     if (abs(t-toutput-dtoutput) <= 1.e-14) then
        stime = MPI_Wtime()
        toutput = t
        noutput = noutput+1
        if (myid == 0) print '(A,i4,A,i8,A,es9.2)', &
             'writing output file ',noutput,', step = ',it,',  t =',t
        call output_mpi(u, t, nxloc, nyloc, myid, pcoords, nx, ny, pdims(1), pdims(2), noutput)
        ftime = MPI_Wtime()
        iotime = iotime + ftime-stime
     end if

  end do


  ! output final solution 
  stime = MPI_Wtime()
  toutput = t
  noutput = noutput+1
  if (myid == 0) print '(A,i4,A,i8,A,es9.2)', &
       'writing output file ',noutput,', step = ',it,',  t =',t
  call output_mpi(u, t, nxloc, nyloc, myid, pcoords, nx, ny, pdims(1), pdims(2), noutput)
  ftime = MPI_Wtime()
  iotime = iotime + ftime-stime
  
  ! output runtime
  if (myid == 0) then
     print *, ' total initialization time = ', inittime
     print *, ' total input/output time   = ', iotime
     print *, ' total communication time  = ', commtime
     print *, ' total simulation time     = ', runtime
  endif

  ! free memory
  call free_exchange()
  deallocate(u, v1, v2, v3)
  deallocate(v2E, v3N, v1W, v1S)

  call MPI_Finalize(ierr)

  ! end program
end program advection_mpi
!=================================================================
