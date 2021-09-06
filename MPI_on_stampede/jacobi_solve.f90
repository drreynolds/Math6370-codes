! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6370
! 23 March 2009
!=================================================================


subroutine jacobi_solve(a, b, c, u, r, res, local_N, global_N, &
     delta, maxiter, iters, comm, ierr)
  !-----------------------------------------------------------------
  ! Description: 
  !   Solves the linear system 
  !           a(k)*u(k-1) + b(k)*u(k) + c(k)*u(k+1) = r(k),
  !   using a parallelized Jacobi iterative solver.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none

  integer, intent(in) :: local_N, global_N, comm, maxiter
  integer, intent(out) :: ierr
  double precision, intent(in), dimension(local_N) :: a, b, c, r
  double precision, dimension(local_N) :: u, res
  double precision, intent(in) :: delta
  integer, intent(out) :: iters
  double precision :: resid2
  integer :: k

  !======= Internals ============
  
  ! compute initial linear residual
  call linresid(a, b, c, u, r, res, resid2, local_N, global_N, comm, ierr)
  if (ierr /= 0) then
     print *,' jacobi_solve error: linresid failed!'
     return
  end if

  ! iterate until resid2 < delta
  do iters=0,maxiter
     
     ! check for convergence 
     if (resid2 < delta)  exit

     ! update u = u - diag(b)^{-1}*res
     u = u - res/b
     
     ! compute linear residual
     call linresid(a, b, c, u, r, res, resid2, local_N, global_N, comm, ierr)
     if (ierr /= 0) then
        print *,' jacobi_solve error: linresid failed!'
        return
     end if
  end do

  if (iters == maxiter) then
     print *,'  jacobi_solve warning: reached maximum iteration limit!'
  end if

  return

end subroutine jacobi_solve
!=================================================================
