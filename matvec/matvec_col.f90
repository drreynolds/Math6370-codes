! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 7 January 2009
!=================================================================


program MatVec
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the product of an m*n matrix and an n-vector.
  !-----------------------------------------------------------------
  !======= Inclusions ===========

  !======= Declarations =========
  implicit none
  integer :: m, n, i, j
  double precision, allocatable :: A(:,:), x(:), b(:)
  double precision :: norm2, stime, ftime

  !======= Internals ============

  ! input the size of the system
  write(*,*) 'We will multiply an m*n matrix by an n-vector'
  write(*,*) '   enter m'
  read(*,*) m
  write(*,*) '   enter n'
  read(*,*) n
  if ((m < 1) .or. (n < 1)) then
     write(*,*)  ' illegal input, m =',m,' and n =',n,&
          ' must both be >= 1'
     stop
  endif
  
  ! initialize the matrix and vectors
  allocate(A(m,n),x(n),b(m))
  do j=1,n
     do i=1,m
        A(i,j) = 1.d0/(1.d0 + (i - j)**2)
     enddo
     x(j) = 1.d0
  enddo

  ! start timer
  call get_time(stime)

  ! compute matrix-vector product (column-based version)
  b = 0.d0
  do j=1,n
     do i=1,m
        b(i) = b(i) + A(i,j)*x(j)
     enddo
  enddo

  ! stop timer
  call get_time(ftime)

  ! output 2-norm of product and runtime to screen
  norm2 = 0.d0
  do i=1,m,1
     norm2 = norm2 + b(i)**2
  enddo
  write(*,*) '       matrix size =',m,' x ',n
  write(*,*) ' 2-norm of product =',sqrt(norm2)
  write(*,*) '           runtime =',ftime-stime

  ! output product to file
  open(101,file='b.txt',form='formatted')
  do i=1,m,1
     write(101,'(es22.15)') b(i)
  enddo
  close(101)
  
  ! free vectors
  deallocate(A,x,b)

end program MatVec
!=================================================================
