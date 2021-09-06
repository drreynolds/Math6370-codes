! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


subroutine output_mpi(u, t, nxloc, nyloc, pid, pcoords, &
                      nx, ny, px, py, noutput)
  !===============================================================
  ! Description: 
  !    Writes current solution to disk in parallel, with a 
  !    separate file for each MPI process.
  !===============================================================
  ! inclusions
  implicit none

  ! declarations
  integer, intent(in) :: nxloc, nyloc, pid, pcoords(2)
  integer, intent(in) :: nx, ny, px, py, noutput
  real*8,  intent(in) :: t, u(nxloc,nyloc)
  integer :: i, j
  character*50 :: outname
  
  ! internals
  
  ! set output file name
  ! Note: we reserve the first set of digits for the MPI process (unused here)
  write(outname,'(5Hu_sol,f4.3,f4.3)') float(pid)/1000.0 , noutput/1000.0

  ! open file
  open(102,file=outname)

  ! write data set parameters for this process
  write(102,*) nxloc
  write(102,*) nyloc
  write(102,*) pcoords(1)
  write(102,*) pcoords(2)
  write(102,*) t

  ! output the solution values and close the data set
  do j=1,nyloc
     do i=1,nxloc
        write(102,*) u(i,j)
     end do
  end do
  close(102)

  ! now output a metadata file, containing general run information
  if (pid == 0) then
     open(103,file='u_sol.txt')
     write(103,*) nx
     write(103,*) ny
     write(103,*) px
     write(103,*) py
     write(103,*) noutput
     close(103)
  end if

  return

  ! end program
end subroutine output_mpi
!=================================================================
