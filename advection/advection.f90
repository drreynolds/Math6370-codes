! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!-----------------------------------------------------------------
! Description: 
!    Evolves the first-order 2D wave equations in time.
!=================================================================


program advection
  !===============================================================
  ! inclusions
  implicit none

  ! interfaces
  interface
     subroutine output(u,t,nx,ny,noutput)
       integer, intent(in) :: nx, ny, noutput
       real*8,  intent(in) :: t, u(nx,ny)
     end subroutine output
     subroutine initialize(u,v1,v2,v3,c,dx,dy,nx,ny)
       integer, intent(in)  :: nx, ny
       real*8,  intent(in)  :: c, dx, dy
       real*8, dimension(nx,ny), intent(out) :: u, v1, v2, v3
     end subroutine initialize
  end interface

  ! declarations
  integer :: nx, ny, nt, i, j, noutput, it
  real*8 :: tstop, c, dtoutput, dx, dy, t, toutput, dt
  real*8 :: stime, ftime, v2_E, v2_W, v3_N, v3_S, v1_W, v1_E, v1_S, v1_N
  real*8 :: runtime, iotime, inittime
  real*8, dimension(:,:), allocatable :: u, v1, v2, v3
  namelist /inputs/ nx, ny, nt, tstop, c, dtoutput
  
  ! internals
  
  call get_time(stime)
  ! read problem parameters from input file (should be in this order):
  !    nx - number of nodes in x-direction
  !    ny - number of nodes in y-direction
  !    nt - number of time steps
  !    tstop - final time (will stop at nt or stop, whichever is 1st)
  !    c - wave speed
  !    dtoutput - time frequency for outputting solutions
  open(100,file='input.txt')
  read(100,inputs)
  close(100)
  print *, '  '
  print *, 'Running wave problem:'
  print *, '   nx =',nx,',  ny =',ny
  print *, '   nt =',nt,',  tstop =',tstop
  print *, '   c =',c
  print *, '   dtoutput =',dtoutput
  print *, '  '
  call get_time(ftime)
  iotime = ftime-stime

  ! allocate arrays
  call get_time(stime)
  allocate(u(nx,ny),v1(nx,ny),v2(nx,ny),v3(nx,ny))

  ! set grid spacing
  dx = 1.d0/nx
  dy = 1.d0/ny

  ! set initial conditions
  call initialize(u,v1,v2,v3,c,dx,dy,nx,ny)
  call get_time(ftime)
  inittime = ftime-stime

  ! set initial time, output initial solution
  t = 0.d0
  toutput = 0.d0
  noutput = 0
  print '(A,i4,A,i8,A,es9.2)', 'writing output file ',noutput, &
       ', step = ',0,',  t =',t
  call get_time(stime)
  call output(u,t,nx,ny,noutput)
  call get_time(ftime)
  iotime = iotime+ftime-stime

  ! set time step
  dt = min(dx/c/5.d1, dy/c/5.d1)

  ! start time stepping
  runtime = 0.d0
  do it=1,nt

     ! start timer
     call get_time(stime)

     ! first update v1 to get to half time step
     !    get relevant values for this location
     do j=1,ny
        do i=1,nx

           ! access relevant components of v2 and v3
           if (i==nx) then
              v2_E = v2(1,j)
           else
              v2_E = v2(i+1,j)
           end if
           v2_W = v2(i,j)
           if (j==ny) then
              v3_N = v3(i,1)
           else
              v3_N = v3(i,j+1)
           end if
           v3_S = v3(i,j)

           ! update v1
           v1(i,j) = v1(i,j) + c*dt/dx*(v2_E - v2_W) + c*dt/dy*(v3_N - v3_S)

        end do
     end do

     ! next update v2 & v3 to get to full time step
     !    get relevant values for this location
     do j=1,ny
        do i=1,nx

           ! access relevant components of v1
           if (i==1) then
              v1_W = v1(nx,j)
           else
              v1_W = v1(i-1,j)
           end if
           v1_E = v1(i,j)
           if (j==1) then
              v1_S = v1(i,ny)
           else
              v1_S = v1(i,j-1)
           end if
           v1_N = v1(i,j)

           ! update v2 and v3
           v2(i,j) = v2(i,j) + c*dt/dx*(v1_E - v1_W)
           v3(i,j) = v3(i,j) + c*dt/dy*(v1_N - v1_S)

        end do
     end do

     ! update and plot solution
     u = u + dt*v1

     ! update time
     t = t + dt
     call get_time(ftime)
     runtime = runtime+ftime-stime

     ! stop simulation if we've reached tstop
     if (t >= tstop)  exit

     ! output solution periodically
     if (abs(t-toutput-dtoutput) <= 1.e-14) then
        call get_time(stime)
        toutput = t
        noutput = noutput+1
        print '(A,i4,A,i8,A,es9.2)', 'writing output file ',noutput, &
             ', step = ',it,',  t =',t
        call output(u,t,nx,ny,noutput)
        call get_time(ftime)
        iotime = iotime+ftime-stime
     end if

  end do


  ! output final solution 
  call get_time(stime)
  toutput = t
  noutput = noutput+1
  print '(A,i4,A,i8,A,es9.2)', 'writing output file ',noutput, &
       ', step = ',it,',  t =',t
  call output(u,t,nx,ny,noutput)
  call get_time(ftime)
  iotime = iotime+ftime-stime

  ! output runtime
  print *, ' total initialization time = ', inittime
  print *, ' total input/output time   = ', iotime
  print *, ' total simulation time     = ', runtime

  ! free memory
  deallocate(u,v1,v2,v3)

  ! end program
end program advection
!=================================================================
