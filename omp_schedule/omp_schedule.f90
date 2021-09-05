! -*- Mode: Fortran90; -*-
!
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370/6370
! 7 February 2015
!=================================================================


program OMPSchedule_Example
  !-----------------------------------------------------------------
  ! Description: 
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer, parameter :: n=23
  integer :: i
  integer, external :: omp_get_thread_num

  !======= Internals ============
  
  !$omp parallel do schedule(runtime)
  do i=1,n
     print *, 'iteration ',i,' performed by thread',omp_get_thread_num()
  enddo
  !$omp end parallel do
  

end program OMPSchedule_Example
!=================================================================
