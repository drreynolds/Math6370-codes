! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================


subroutine output(u,t,nx,ny,noutput)
  !===============================================================
  ! Description:
  !    Writes current solution to disk.
  !===============================================================
  ! inclusions
  implicit none

  ! declarations
  integer, intent(in) :: nx, ny, noutput
  real*8,  intent(in) :: t, u(nx,ny)
  integer :: i, j
  character*50 :: outname

  ! internals

  ! set output file name
  ! Note: we reserve the first set of digits for the MPI process (unused here)
  write(outname,'(5Hu_sol,f4.3,f4.3)') 0.0, noutput/1000.0

  ! open file
  open(102,file=outname)

  ! write data set parameters
  ! Note: the two 0's will be used for the MPI process location (unused here)
  write(102,*) nx
  write(102,*) ny
  write(102,*) 0
  write(102,*) 0
  write(102,*) t

  ! output the solution values and close the data set
  do j=1,ny
     do i=1,nx
        write(102,*) u(i,j)
     end do
  end do
  close(102)


  ! now output a metadata file, containing general run information
  ! Note: the two 1's will be used for the MPI process dimensions (unused here)
  open(103,file='u_sol.txt')
  write(103,*) nx
  write(103,*) ny
  write(103,*) 1
  write(103,*) 1
  write(103,*) noutput
  close(103)

  return

  ! end program
end subroutine output
!=================================================================
