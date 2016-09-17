      subroutine laplace_approximation(n, c, K, y, sig, D)
!
      implicit none
      integer, intent(in)           :: n
      double precision, intent(in)  :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)
! 
      double precision, allocatable :: yold(:), S(:, :), DIAG(:, :), CURR(:, :), work(:)
      integer,          allocatable :: ipiv(:)
!
      double precision, parameter :: thresh = 0.0000001
      integer, parameter :: niter = 100000
      integer :: i, info, iter = 1
!
      external DGETRF
      external DGETRI

      allocate(yold(n), S(n, n), DIAG(n, n), CURR(n, n), work(n))
      allocate(ipiv(n))
!
      DIAG = 0
      D = 0
      forall(i = 1:n) DIAG(i, i) = 1
      yold = 1.0
!
      do while (sum(abs(y - yold)) > thresh .and. iter < niter)
          yold = y
          forall(i = 1:n) sig(i) = 1 / (1 + exp(y(i)))
          forall(i = 1:n) D(i, i) = y(i) * (1 - y(i))
          CURR = DIAG + MATMUL(D, K)
          call DGETRF(n, n, CURR, n, ipiv, info)
          call DGETRI(n, CURR, n, ipiv, work, n, info)
          y = MATMUL(MATMUL(K, CURR), MATMUL(D, y) + c - sig)
          iter = iter + 1
      end do

      deallocate(yold, S, DIAG, CURR, work, ipiv)
      return
      end subroutine

