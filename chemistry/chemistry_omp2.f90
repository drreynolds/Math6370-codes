! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 6370
! 27 March 2011
!=================================================================


double precision function TIMER()
  !-----------------------------------------------------------------
  ! Function to use the best available timer
  !-----------------------------------------------------------------
  double precision :: dtime
  double precision, external :: omp_get_wtime
  call get_time(TIMER)
  !$TIMER = omp_get_wtime()
end function TIMER


program ChemicalEquilibrium_OpenMP
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
     subroutine chem_solver_orphan(n,T,u,v,w,lam,eps,maxit,mycount,err)
       integer, intent(in) :: n, maxit
       integer, intent(out) :: mycount, err
       double precision, intent(in) :: T(n), lam, eps
       double precision, intent(inout) :: u(n), v(n), w(n)
     end subroutine chem_solver_orphan
  end interface

  !======= Declarations =========
  integer :: n, maxit, mycount, numprocs, err=0, locerr
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: lam, eps, stime, ftime
  integer, external :: omp_get_num_threads, omp_get_thread_num
  double precision, external :: TIMER

  !======= Internals ============
  
  ! set solver input parameters
  maxit = 10000000
  lam = 1.0d-2
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
  stime = TIMER()

  ! start OpenMP parallelism
  !$omp parallel default(shared) private(mycount,locerr)

  ! output parallelism information
  !$omp master
  !$ numprocs = omp_get_num_threads()
  !$ write(*,'(A,i5,A)') ' starting OpenMP with ', numprocs,' processes'
  !$omp end master

  ! call solver 
  call chem_solver_orphan(n,T,u,v,w,lam,eps,maxit,mycount,locerr)

  ! accumulate errors
  !$omp critical
  err = err + locerr
  !$omp end critical

  ! output loop partitioning information
  !#omp critical
  !$ write(*,'(2(A,i5),A)') ' thread ',omp_get_thread_num(),' computed ',mycount,' iterations'
  !#omp end critical

  !$omp end parallel
  if (err > 0) then
     write(*,*) '*** WARNING ***'
     write(*,*) '  See above messages for errors in solver'
  endif

  ! stop timer
  ftime = TIMER()

  ! output solution time
  write(*,*) '     runtime =',ftime-stime


  ! free temperature and solution arrays
  deallocate(T,u,v,w)

end program ChemicalEquilibrium_OpenMP
!=================================================================
