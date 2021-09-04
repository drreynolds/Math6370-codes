! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


subroutine initialize(u,v1,v2,v3,c,dx,dy,nx,ny)
  !===============================================================
  ! Description: 
  !    Sets the initial conditions into u, v1, v2, v3.
  !===============================================================
  ! inclusions
  implicit none

  ! declarations
  integer, intent(in)  :: nx, ny
  real*8,  intent(in)  :: c, dx, dy
  real*8, dimension(nx,ny), intent(out) :: u, v1, v2, v3
  real*8 :: xspan_c(nx), xspan_h(nx), yspan_c(ny), yspan_h(ny)
  real*8 :: x_c, x_h, y_c, y_h
  integer :: i, j
  
  ! internals
  
  ! set mesh points
  do i=1,nx
     xspan_c(i) = dx/2.d0 + (i-1)*dx
     xspan_h(i) = (i-1)*dx
  end do

  do j=1,ny
     yspan_c(j) = dy/2.d0 + (j-1)*dy
     yspan_h(j) = (j-1)*dy
  end do


  ! set initial condition for solution and derivatives
  do j=1,ny

     ! y locations
     y_c = yspan_c(j)
     y_h = yspan_h(j)

     do i=1,nx

        ! x locations
        x_c = xspan_c(i)
        x_h = xspan_h(i)

        ! set initial conditions on u_x, u_y [gaussian blob]
        !    u(0,x,y) = exp(-100*((x-1/3)^2+(y-1/2)^2))
        !    c*u_x(0,x,y) = -200*c*(x-1/3)*exp(-100*((x-1/3)^2+(y-1/2)^2))
        !    c*u_y(0,x,y) = -200*c*(y-1/2)*exp(-100*((x-1/3)^2+(y-1/2)^2))
        u(i,j) = exp(-1.d2*((x_c-1.d0/3.d0)**2+(y_c-0.5d0)**2))
        v1(i,j) = 0.d0
        v2(i,j) = -2.d2*c*(x_h-1.d0/3.d0)*exp(-1.d2*((x_h-1.d0/3.d0)**2+(y_c-0.5d0)**2))
        v3(i,j) = -2.d2*c*(y_h-0.5d0)*exp(-1.d2*((x_c-1.d0/3.d0)**2+(y_h-0.5d0)**2))

     end do
  end do
  
  return

  ! end program
end subroutine initialize
!=================================================================
