! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


program Global_Min
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes a global minimum to the function
  !        f(x,y) = exp(sin(50x)) + sin(60exp(y)) + sin(70sin(x))
  !                + sin(sin(80y)) - sin(10(x+y)) + (x^2+y^2)/4
  !    We start with a simple search algorithm that quickly finds 
  !    100 suitable starting points for local minimization 
  !    algorithms (steepest descent).  Each of these starting
  !    points are then examined thoroughly to find the nearest 
  !    local minimum.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer, parameter :: nx = 1000
  integer, parameter :: ny = 1000
  integer, parameter :: npts = 100
  integer, parameter :: maxits = 100000000
  integer :: i, ix, iy, k, np, idx(1)
  double precision :: f, fx, fy, x, y, dx, dy, stime, ftime
  double precision :: pt(2), tstpt(2), df(2), bestpt(2)
  double precision :: cutoff, curval, fval, gamma, bestval
  double precision :: searchpts(2,npts), searchvals(npts)

  !======= Internals ============
  
  ! set the integrand function
  f(x,y) = exp(sin(5.d1*x)) + sin(6.d1*exp(y)) + sin(7.d1*sin(x)) &
         + sin(sin(8.d1*y)) - sin(1.d1*(x+y)) + (x**2+y**2)/4.d0

  ! set the jacobian functions
  fx(x,y) = exp(sin(5.d1*x))*cos(5.d1*x)*5.d1 &
          + cos(7.d1*sin(x))*7.d1*cos(x)      &
          - cos(1.d1*(x+y))*1.d1 + x/2.d0
  fy(x,y) = cos(6.d1*exp(y))*6.d1*exp(y)      &
          + cos(sin(8.d1*y))*cos(8.d1*y)*8.d1 &
          - cos(1.d1*(x+y))*1.d1 + y/2.d0


  ! start timer
  call get_time(stime)

  ! set subinterval widths
  ! (note: we know that the minimum is inside the box [-5,5]^2)
  dx = 10.d0/(nx-1)
  dy = 10.d0/(ny-1)

  print *, 'initial search over mesh'
  ! check initial mesh, saving best npts points
  cutoff = 1.d12     ! initialize to very large number
  np = 0             ! no points in queue yet
  do i=1,nx*ny

     ! set mesh location
     iy = (i-1)/nx
     ix = i - iy*nx
     pt(1) = -5.d0 + (ix-1)*dx
     pt(2) = -5.d0 + iy*dy

     ! evaluate current point
     curval = f(pt(1),pt(2))

     ! if current value is below cutoff, add to list of points
     if (curval < cutoff) then
        
        ! if list has room, just add this point
        if (np < npts-1) then
           np = np+1                ! increment counter
           searchpts(:,np) = pt     ! add point to list
           searchvals(np) = curval  ! add value at point
           
        ! if this is the last empty slot in the list, add point 
        ! and set cutoff
        else if (np == npts-1) then
           np = np+1                    ! increment counter
           searchpts(:,np) = pt         ! add point to list
           searchvals(np) = curval      ! add value at point
           cutoff = maxval(searchvals)  ! set cutoff as worst value
           
        ! otherwise replace an inferior entry and update cutoff
        else
           idx = maxloc(searchvals)     ! index of worst pt in list
           searchpts(:,idx(1)) = pt     ! replace point
           searchvals(idx(1)) = curval  ! replace value
           cutoff = maxval(searchvals)  ! update cutoff

        end if
        
     end if
  end do

  print *, 'performing minimizations over best ',npts,' points'
  ! We have our list of the best npts test points.  We now do 
  ! local minimization (via steepest descent) around each of 
  ! these to better refine
  bestval = 1.d12     ! initialize to very large number
  do i = 1,npts
     
     ! extract the current point and its function value
     pt = searchpts(:,i)
     fval = searchvals(i)

     ! perform a steepest descent minimization at this point
     do ix = 1,maxits

        ! compute gradient of f at this point
        df(1) = fx(pt(1),pt(2))
        df(2) = fy(pt(1),pt(2))

        ! set the initial linesearch step size
        gamma = 1.d0/sqrt(df(1)**2 + df(2)**2)

        ! perform back-tracking line search for gamma
        do k = 1,50

           ! set test point and calculate function value
           tstpt = pt - gamma*df
           curval = f(tstpt(1),tstpt(2))

           ! if test point successful, exit; otherwise reduce gamma
           if (curval < fval) then
              exit
           else 
              gamma = 0.5d0*gamma
           end if

        end do
        
        ! check for stagnation/convergence
        if (maxval(abs((pt-tstpt))) < 1e-13)  exit

        ! update point with current iterate
        pt = tstpt
        fval = curval

     end do

     ! if current value is better than "best" so far, update best
     if (fval < bestval) then
        bestpt = pt
        bestval = fval
        print *,' new best-guess has value ',bestval
     end if

  end do

  ! stop timer
  call get_time(ftime)

  ! output computed minimum and corresponding point
  write(*,*) ' computed minimum =',bestval
  write(*,*) '            point =',bestpt(1),bestpt(2)
  write(*,*) '          runtime =',ftime-stime
  

end program Global_Min
!=================================================================
