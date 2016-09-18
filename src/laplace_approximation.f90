
      subroutine newton_step(n, c, K, y, sig, D, DIAG)
! compute Newton update
! references: Barber p406, Equation~19.5.19; Rasmussen p43, Equation~3.18
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
! compute the Laplace approximation around the mode
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
      double precision, allocatable :: CURR(:, :), work(:)
      integer,          allocatable :: ipiv(:)
      integer :: info

      external DGETRF
      external DGETRI

      allocate(yold(n), DIAG(n, n))
      allocate(CURR(n, n), work(n))
      allocate(ipiv(n))
!
      DIAG = 0
      D = 0
      forall(i = 1:n) DIAG(i, i) = 1
      yold = 1.0
      y = 0.0
      !
      do while (sum(abs(y - yold)) > thresh .and. iter < niter)
          yold = y
! compute class mapping for y        
          forall(i = 1:n) sig(i) = 1 / (1 + exp(-y(i)))
! compute Hessian of log-likelihood (ONLY IN THIS CASE WHERE THE LOG-TRANSFORM IS USED!)          
          forall(i = 1:n) D(i, i) = sig(i) * (1 - sig(i))
          CURR = MATMUL(D, K) + DIAG
          call DGETRF(n, n, CURR, n, ipiv, info)
! invert CURR          
          call DGETRI(n, CURR, n, ipiv, work, n, info)
! calculate Newton step
          y = MATMUL(MATMUL(K, CURR), MATMUL(D, y) + c - sig)
          ! call newton_step(n, c, K, y, sig, D, DIAG)
          iter = iter + 1
      end do
!
      deallocate(yold, DIAG)
      deallocate(CURR, work, ipiv)
      return
      end subroutine
