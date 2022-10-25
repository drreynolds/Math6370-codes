! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================


subroutine get_time(tval)
  !-----------------------------------------------------------------
  ! Description:
  !    Calls the date_and_time routine to output the current time
  !    of day in seconds, and returns with a negative value if the
  !    system_clock routine does not provide complete information.
  !
  ! Arguments:
  !    tval - double (output), current time in seconds
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none
  integer :: c, r, m
  double precision :: tval

  !======= Internals ============

  call system_clock(COUNT = c, COUNT_RATE = r, COUNT_MAX = m)

  ! check that system_clock gives complete information
  if (c < 0) then
     tval = -1
     return
  endif

  ! combine outputs into the result
  tval = 1.d0*c/r

  return

end subroutine get_time
!=================================================================
