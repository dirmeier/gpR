      subroutine laplace_approximation(n, c, K, y, sig, D)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)i

      double precision :: yold(n), sig(n), S(n, n), DIAG(n, n)
      double precision, parameter :: thresh = 0.0000001
      integer, parameter niter :: 100000
      integer :: i, iter = 1

      external DGETRF
      external DGETRI

      DIAG = 0
      forall(i = 1:n) DIAG(i, i) = 1

      yold = 1.0
      do while (sum(abs(y - yold)) > thresh && iter < niter)
          sig = sigmoid(y, n)
          D = hessian(sig, n)
          yold = y
          K = DIAG + MATMUL(D, K)
          call DGETRF(n, n, Ainv, n, ipiv, info)
          call DGETRI(n, Ainv, n, ipiv, work, n, info)
          y = MATMUL(K, )
          iter = iter + 1
      end do

      return
      end subroutine

       while (base::mean(base::abs(y - y.old)) > 10e-5 && iter <= niter)
    {
      # compute class mapping for y
      sig <- as.vector(sigmoid(y))
      # compute Hessian of log-likelihood (ONLY IN THIS CASE WHERE THE LOG-TRANSFORM IS USED!)
      D <- diag(sig*(1-sig))
      y.old <- y
      # compute Newton update
      #  (see references: Barber p406, Equation~19.5.19; Rasmussen p43, Equation~3.18)
      y <- base::as.vector(K.train %*%
                             base::solve(diagn + D %*% K.train) %*%
                             (D%*%y + c.train  - sig))
      iter <- iter + 1
    }

      pure function hessian(y, n)
      integer, intent(in) :: n
      double precision, intent(in) :: y(n)
      double precision, intent(out) :: hessian(n, n)
      integer :: i
      hessian (:, :) = 0.0
      do i = 1, n
          hessian(i, i) = y * (1 - y)
      end do
      end function

      pure function sigmoid(y, n)
      integer intent(in):: n
      double precision, intent(in) ::  y(n)
      double precision, intent(out) :: sigmoid(n)
      integer :: i
      do i = 1, n
          sigmoid(i) = 1 / (1 + exp(y(i)))
      end do
      end function
