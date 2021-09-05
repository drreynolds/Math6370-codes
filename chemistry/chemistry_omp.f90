! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 7 January 2009
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
     subroutine chem_solver(T,u,v,w,lam,eps,maxit,its,res)
       integer, intent(in) :: maxit
       integer, intent(out) :: its
       double precision, intent(in) :: T, lam, eps
       double precision, intent(out) :: res
       double precision, intent(inout) :: u, v, w
     end subroutine chem_solver
  end interface

  !======= Declarations =========
  integer :: n, i, maxit, its, numprocs, mycount, err=0
  double precision, allocatable :: T(:), u(:), v(:), w(:)
  double precision :: lam, eps, res, stime, ftime
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
  !$omp parallel default(shared) private(mycount,its,res)

  ! output parallelism information
  !$omp single
  !$ numprocs = omp_get_num_threads()
  !$ write(*,'(A,i5,A)') ' starting OpenMP with ', numprocs,' processes'
  !$omp end single

  ! call solver over n intervals
  mycount = 0
  !$omp do schedule(static, n/numprocs) reduction(+:err)
  do i=1,n,1
     call chem_solver(T(i),u(i),v(i),w(i),lam,eps,maxit,its,res)
     if (res < eps) then
!!$        write(*,'(2(A,i8))') '    i =',i,',  its =',its
     else
        write(*,'(2(A,i8),4(A,es9.2))') '    error: i =',i,',  its =',its, &
             ',  res =',res,',  u =',u(i),',  v =',v(i),',  w =',w(i)
        err = err + 1
     endif

     ! update thread increment
     mycount = mycount + 1

  enddo
  !$omp end do

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
