c=================================================================
c     Daniel R. Reynolds
c     SMU Mathematics
c     Math 4370 / 6370
c=================================================================


      program DotProd
c-----------------------------------------------------------------
c     Example routine to compute the dot-product of two vectors.
c-----------------------------------------------------------------
c======= Inclusions ===========

c======= Declarations =========
      implicit none

      integer, parameter :: MAX_INTEGER_DIGITS = 10
      integer :: i, n
      character(MAX_INTEGER_DIGITS) :: n_string
      double precision, allocatable :: a(:), b(:)
      double precision :: sum, stime, ftime
      double precision :: alloctime, inittime, runtime

c======= Internals ============

c     get n from the command line
      call getarg(1, n_string)

c     ensure that an argument was passed in
      if ( n_string == '' ) then
         stop 'Error: function requires one argument (vector length)'
      endif

c     convert n_string to integer, and ensure it's positive
      read (n_string, *) n
      if ( n < 1 ) then
         stop 'Error: vector length must be greater than 0'
      endif

c     allocate the vectors
      call get_time(stime)
      allocate(a(n),b(n))
      call get_time(ftime)
      alloctime = ftime - stime

c     initialize the vectors
      call get_time(stime)
      do i=1,n
         a(i) = 1.d-3*(i)/n
         b(i) = 1.d-3*(n-i)/n
      enddo
      call get_time(ftime)
      inittime = ftime - stime

c     compute dot-product
      call get_time(stime)
      sum = 0.d0
      do i=1,n,1
         sum = sum + a(i)*b(i)
      enddo
      call get_time(ftime)
      runtime = ftime-stime

c     output computed value and runtime
      print *, ' vector length =',n
      print *, '   dot-product =',sum
      print *, '    alloc time =',alloctime
      print *, '     init time =',inittime
      print *, '      run time =',runtime

c     free vectors
      deallocate(a,b)

      end program DotProd
c=================================================================
