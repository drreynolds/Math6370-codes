! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 7 January 2009
!=================================================================


program MatVec_OpenMP
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the product of an m*n matrix and an n-vector.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none

  integer :: m, n, i, j, numprocs
  double precision, allocatable :: A(:,:), x(:), b(:)
  double precision :: norm2, stime, ftime
  integer, external :: omp_get_num_threads
  double precision, external :: omp_get_wtime

  !======= Internals ============
  
  ! input the size of the system
  write(*,*) 'We will multiply an m*n matrix by an n-vector'
  write(*,*) '   enter m'
  read(*,*) m
  write(*,*) '   enter n'
  read(*,*) n
  if ((m < 1) .or. (n < 1)) then
     write(*,*)  ' illegal input, m =',m,' and n =',n,&
          ' must both be positive'
     stop
  endif
  
  ! allocate the matrix and vectors
  allocate(A(m,n),x(n),b(m))

  ! initialize the matrix and vectors
  do j=1,n
     do i=1,m
        A(i,j) = 1.d0/(1+(i-j)**2)
     enddo
  enddo
  b = 0.d0
  x = 1.d0

  ! start timer
  stime = omp_get_wtime()

  ! start OpenMP parallelism
  !$omp parallel default(shared) private(i,j)

  ! output parallelism information, start timer
  !$omp single
  write(*,'(A,i5,A)') ' starting OpenMP with ', omp_get_num_threads(),' processes'
  !$omp end single


  ! compute matrix-vector product
  do j=1,n
     !$omp do
     do i=1,m
        b(i) = b(i) + A(i,j)*x(j)
     enddo
     !$omp end do
  enddo

  ! end OpenMP parallel region
  !$omp end parallel

  ! stop timer
  ftime = omp_get_wtime()

  ! output 2-norm of product and runtime to screen
  norm2 = sum(b*b)
  write(*,*) '       matrix size =',m,' x ',n
  write(*,*) ' 2-norm of product =',sqrt(norm2)
  write(*,*) '           runtime =',ftime-stime

  ! output product to file
  open(101,file='b_omp.txt',form='formatted')
  do i=1,m,1
     write(101,'(es22.15)') b(i)
  enddo
  close(101)
  
  ! free vectors
  deallocate(A,x,b)

end program MatVec_OpenMP
!=================================================================
