     subroutine y_posterior(n, c, K, y, sig, D, DIAG)
!
      implicit none
      integer, intent(in)           :: n
      double precision, intent(in)  :: c(n), K(n, n), sig(n), D(n, n), DIAG(n, n)
      double precision, intent(out) :: y(n)
!
      double precision, allocatable :: CURR(:, :), work(:), bla(:)
      integer,          allocatable :: ipiv(:)
      integer :: info
!
      external DGETRF
      external DGETRI
!
      allocate(CURR(n, n), work(n), bla(N))
      allocate(ipiv(n))
!
      CURR = MATMUL(D, K) + DIAG
      call DGETRF(n, n, CURR, n, ipiv, info)
      call DGETRI(n, CURR, n, ipiv, work, n, info)
      y = MATMUL(MATMUL(K, CURR), MATMUL(D, y) + c - sig)
!
      deallocate(CURR, work, ipiv, bla)
      return
      end subroutine
!
!
      subroutine laplace_approximation(n, c, K, y, sig, D)
!
      implicit none
      integer, intent(in)           :: n
      double precision, intent(in)  :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)
!
      double precision, allocatable :: yold(:), DIAG(:, :)
!
      double precision, parameter :: thresh = 0.00001
      integer, parameter :: niter = 10000
      integer :: i, iter = 1
!
      allocate(yold(n), DIAG(n, n))
!
      DIAG = 0
      D = 0
      forall(i = 1:n) DIAG(i, i) = 1
      yold = 1.0
!
      do while (sum(abs(y - yold)) > thresh .and. iter < niter)
          yold = y
          forall(i = 1:n) sig(i) = 1 / (1 + exp(y(i)))
          forall(i = 1:n) D(i, i) = sig(i) * (1 - sig(i))
          call y_posterior(n, c, K, y, sig, D, DIAG)
          iter = iter + 1
      end do
!
      deallocate(yold, DIAG)
      return
      end subroutine
