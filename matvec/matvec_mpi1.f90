! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU Mathematics
! Math 4370 / 6370
!=================================================================


program MatVec_MPI
  !-----------------------------------------------------------------
  ! Description:
  !    Computes the product of an m*n matrix and an n-vector.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none
  integer :: m, n, i, j, is, ie, numprocs, myid, ierr
  double precision, allocatable :: A(:,:), x(:), myb(:), b(:)
  double precision :: norm2, stime, ftime

  !======= Internals ============

  ! intialize MPI
  call MPI_Init(ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Init =',ierr
     stop
  endif
  call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_size =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Comm_rank =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! input the size of the system
  if (myid == 0) then
     print *, 'We will multiply an m*n matrix by an n-vector'
     print *, '   enter m'
     call flush()
     read(*,*) m
     print *, '   enter n'
     call flush()
     read(*,*) n
  endif

  ! root node sends m and n out to other processors
  call MPI_Bcast(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif
  call MPI_Bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Bcast =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  if ((m < 1) .or. (n < 1)) then
     if (myid == 0) then
        write(0,*)  ' illegal input, m =',m,' and n =',n,&
             ' must both be >= 1'
     endif
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! root node outputs parallelism information to screen
  if (myid == 0) then
     print '(A,i5,A)', ' starting MPI with ', numprocs,' processes'
  endif

  ! allocate the matrix and vectors
  ! (store entire system on every proc -- wasteful)
  allocate(A(m,n),x(n),myb(m),b(m))

  ! initialize the matrix and vectors
  do j=1,n
     do i=1,m
        A(i,j) = 1.d0/(1+(i-j)**2)
     enddo
  enddo
  b = 0.d0
  x = 1.d0

  ! start timer
  stime = MPI_Wtime()

  ! determine this processor's interval
  is = floor(1.d0*m/numprocs)*myid + 1
  ie = floor(1.d0*m/numprocs)*(myid+1)
  if (myid == numprocs-1) then
     ie = m
  endif

  ! compute matrix-vector product
  myb = 0.d0
  do j=1,n
     do i=is,ie
        myb(i) = myb(i) + A(i,j)*x(j)
     enddo
  enddo

  ! root node collects result
  call MPI_Reduce(myb, b, m, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  if (ierr /= 0) then
     write(0,*) ' error in MPI_Reduce =',ierr
     call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  endif

  ! stop timer
  ftime = MPI_Wtime();

  ! output 2-norm of product and runtime to screen
  if (myid == 0) then
     norm2 = 0.d0
     do i=1,m
        norm2 = norm2 + b(i)**2
     enddo
     print *, '       matrix size =',m,' x ',n
     print *, ' 2-norm of product =',sqrt(norm2)
     print *, '           runtime =',ftime-stime

     ! output product to file
     open(101,file='b_mpi1.txt',form='formatted')
     do i=1,m
        write(101,'(es22.15)') b(i)
     enddo
     close(101)
  endif

  ! free vectors
  deallocate(A,x,b,myb)

  ! finalize MPI
  call MPI_Finalize(ierr)

end program MatVec_MPI
!=================================================================
