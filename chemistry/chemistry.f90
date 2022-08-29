! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 4370 / 6370
!=================================================================


program ChemicalEquilibrium
  !-----------------------------------------------------------------
  ! Description:
  !    Computes the equilibrium chemical densities at a number of
  !    spatial locations, given a (random) background temperature
  !    field.  The chemical rate equations and solution strategy
  !    are in the subroutine chem_solver, which is called at every
  !    spatial location.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  implicit none

  !======= Interfaces ===========
  interface
     subroutine chem_solver(T,u,v,w,lam,eps,maxit,its,res)
       integer, intent(in) :: maxit
       integer, intent(out) :: its
       double precision, intent(in) :: T, lam, eps
       double precision, intent(out) :: res
       double precision, intent(inout) :: u, v, w
     end subroutine chem_solver
  end interface

  !======= Declarations =========
  integer :: n, i, maxit, its
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: lam, eps, res, stime, ftime

  !======= Internals ============

  ! set solver input parameters
  maxit = 1000000
  lam = 1.d-2
  eps = 1.d-10

  ! input the number of intervals
  write(*,*) 'Enter the number of intervals (0 quits):'
  read(*,*) n
  if (n < 1) then
     stop
  endif

  ! allocate temperature and solution arrays
  allocate(T(n),u(n),v(n),w(n))

  ! set random temperature field, initial guesses at chemical densities
  call random_number(T)
  u = 0.35d0
  v = 0.1d0
  w = 0.5d0

  ! start timer
  call get_time(stime)

  ! call solver over n intervals
  do i=1,n,1
     call chem_solver(T(i),u(i),v(i),w(i),lam,eps,maxit,its,res)
     if (res < eps) then
        write(*,'(2(A,i8))') '    i =',i,',  its =',its
     else
        write(*,'(2(A,i8),4(A,es9.2))') '    error: i =',i,',  its =',its, &
             ',  res =',res,',  u =',u(i),',  v =',v(i),',  w =',w(i)
        stop
     endif
  enddo

  ! stop timer
  call get_time(ftime)

  ! output solution time
  write(*,*) '     runtime =',ftime-stime

  ! free temperature and solution arrays
  deallocate(T,u,v,w)

end program ChemicalEquilibrium
!=================================================================
