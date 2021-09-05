! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 13 January 2013
!=================================================================


program DotProd_OpenMP
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the dot-product of two vectors.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer :: i, n
  double precision, allocatable :: a(:), b(:)
  double precision :: sum, alloctime, inittime, runtime
  double precision :: stime, ftime
  double precision, external :: omp_get_wtime
  
  !======= Internals ============
  
  ! get the vector length from the user
  print *, 'input vector length n (must be >0):'
  read(*,*) n

  ! ensure that n is positive
  if (n < 1) then
     print *, 'Error: vector length', n, ' must be greater than 0'
     stop
  endif

  ! alloctate the vectors
  stime = omp_get_wtime()
  allocate(a(n),b(n))
  ftime = omp_get_wtime()
  alloctime = ftime-stime

  ! initialize the vectors
  stime = omp_get_wtime()
  do i=1,n
     a(i) = 1.d-3*i/n
  enddo
  do i=1,n
     b(i) = 1.d-3*(n-i)/n
  enddo
  ftime = omp_get_wtime()
  inittime = ftime-stime

  ! perform parallel dot-product
  stime = omp_get_wtime()
  sum = 0.d0
  !$omp parallel do reduction(+:sum)
  do i=1,n,1
     sum = sum + a(i)*b(i)
  enddo
  !$omp end parallel do
  ftime = omp_get_wtime()
  runtime = ftime-stime

  ! output computed value and runtime
  print *, ' vector length =',n
  print *, '   dot-product =',sum
  print *, '    alloc time =',alloctime
  print *, '     init time =',inittime
  print *, '      run time =',runtime

  ! free vectors
  deallocate(a,b)

end program DotProd_OpenMP
!=================================================================
