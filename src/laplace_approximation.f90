!! Fortran routines for calculation of Laplace approximation and Newton updated.
!!

! compute Newton update
! references: Barber p406, Equation~19.5.19; Rasmussen p43, Equation~3.18
      subroutine newton_step(n, c, K, y, sig, D, DIAG)
      use gpr_util
      implicit none
      integer, intent(in)           :: n
      double precision, intent(in)  :: c(n), K(n, n), sig(n), D(n, n), DIAG(n, n)
      double precision, intent(out) :: y(n)

      double precision, allocatable :: CURR(:, :)

      allocate(CURR(n, n), work(n))
      allocate(ipiv(n))
! calculate the Newton update
      CURR = MATMUL(D, K) + DIAG
      call solve(n, CURR)
      y = MATMUL(MATMUL(K, CURR), MATMUL(D, y) + c - sig)

      deallocate(CURR)
      return
      end subroutine

! compute the Laplace approximation around the mode
      subroutine laplace_approximation(n, c, K, y, sig, D)
      implicit none
      integer, intent(in)           :: n
      double precision, intent(in)  :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)
!
      double precision, allocatable :: yold(:), DIAG(:, :)
!
      double precision, parameter :: thresh = 0.00001
      integer, parameter :: niter = 10000
      integer :: i, iter
!
      allocate(yold(n), DIAG(n, n))
!
      DIAG = 0
      D = 0
      forall(i = 1:n) DIAG(i, i) = 1
      yold = 1.0
      y = 0.0
      iter = 1
! do IRLS until convergence
      do while (sum(abs(y - yold)) > thresh .and. iter < niter)
          yold = y
! compute class mapping for y
          forall(i = 1:n) sig(i) = 1 / (1 + exp(-y(i)))
! compute Hessian of log-likelihood (ONLY IN THIS CASE WHERE THE LOG-TRANSFORM IS USED!)
          forall(i = 1:n) D(i, i) = sig(i) * (1 - sig(i))
! calculate Newton step
          call newton_step(n, c, K, y, sig, D, DIAG)
          iter = iter + 1
      end do
!
      deallocate(yold, DIAG)
      return
      end subroutine

