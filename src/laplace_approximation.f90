      subroutine laplace_approximation(n, c, K, y, sig, D)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)

      double precision :: yold(n), S(n, n), DIAG(n, n), CURR(n, n), work(n)
      double precision, parameter :: thresh = 0.0000001
      integer, parameter :: niter = 100000
      integer :: ipiv(n), i, info, iter = 1

      external DGETRF
      external DGETRI

      DIAG = 0
      forall(i = 1:n) DIAG(i, i) = 1

      yold = 1.0
      do while (sum(abs(y - yold)) > thresh .and. iter < niter)
          sig = sigmoid(y, n)
          D = hessian(sig, n)
          yold = y
          CURR = DIAG + MATMUL(D, K)
          call DGETRF(n, n, CURR, n, ipiv, info)
          call DGETRI(n, CURR, n, ipiv, work, n, info)
          y = MATMUL(MATMUL(K, CURR), MATMUL(D, y) + c - sig)
          iter = iter + 1
      end do

      return
      end subroutine

      function hessian(y, n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: y(n)
      double precision :: hessian(n, n)
      integer :: i
      hessian (:, :) = 0.0
      do i = 1, n
          hessian(i, i) = y * (1 - y)
      end do
      end function

      function sigmoid(y, n)
      implicit none
      integer, intent(in):: n
      double precision, intent(in) :: y(n)
      double precision :: sigmoid(n)
      integer :: i
      do i = 1, n
          sigmoid(i) = 1 / (1 + exp(y(i)))
      end do
      end function
