!! Fortran routines for calculation of binomial class labels
!!

! invert a matrix
      subroutine solve(n, D)
      implicit none
      integer, intent(in)           :: n
      double precision, intent(inout) :: D(n, n)
! worker arrays
      double precision, allocatable :: work(:)
      integer,          allocatable :: ipiv(:)
      integer :: info
! load LAPACK routines for matrix inversion
      external DGETRF
      external DGETRI
! allocate arrays
      allocate(work(n), ipiv(n))
! do matrix inversion
      call DGETRF(n, n, D, n, ipiv, info)
      call DGETRI(n, D, n, ipiv, work, n, info)
! free matrices
      deallocate(work, ipiv)
      return
      end subroutine


! calculate the logistic normal integral using the error function
      pure function lognint(m, var) result(e)
      implicit none
      double precision, intent(in) :: m, var
      double precision, parameter :: PI = 3.1415926535897932
      double precision            :: e, lam
      lam = sqrt(PI / 4)
      e = .5 + .5 * erf(lam * m / sqrt(1.0 + (2 * (lam ** 2.0) * var)))
      end function lognint

! draw from a Gaussian distribution using polar Box-Muller transform
      function rnorm(m, sd)
      implicit none
      double precision, intent(in) :: m, sd
      double precision             :: x1, x2, w, y1, y2

      w = 1.0
      do while( i < 5 )
          x1 = 2.0 * rand() - 1.0;
          x2 = 2.0 * rand() - 1.0;
          w = x1 * x1 + x2 * x2;
      end do

      end function rnorm

! compute the binomial class labels
      subroutine predict_binomial(na, ntrain, nnew, K, ctrain, sigtrain, xnew, trainidx, newidx, me, ps, pm, Kn)
      implicit none
      integer, intent(in)           :: na, ntrain, nnew
      integer, intent(in)           :: trainidx(ntrain), testidx(nnew)
      double precision, intent(in)  :: ct(ntrain), K(na, na), sig(ntrain), xnew(nnew),
      double precision, intent(out) :: me(ntest), ps(nnew), pm(nnew), Kn(nnew, nnew)
!
      double precision :: knn
      double precision, allocatable :: cm(:), Ktt(:, :), D(:, :), Knt(:), m(:), Cu(:, :)
      integer :: i
!
      allocate(cm(ntrain), Ktt(ntrain, ntrain), D(ntrain, ntrain), Knt(n), m(n), Cu(ntrain, ntrain))
!
      Ktt = K(1:trainidx, 1:trainidx)
      cm = ct - sig
      D = 0.0
      forall(i = 1:n) D(i, i) = sig(i) * (1 - sig(i))
      do i = 1, n
          Knt = K(1:trainidx, newidx(i))
          knn = k(newidx(i), newidx(i))
          Cu  = Ktt + D
          solve(n, Cu)
          m   = MATMUL(Knt, cm)
          var = knn - MATMUl(MATMUL(Knt, Cu), Knt)
          me  = lognint(m, var)
          pre =  1 / (1 + exp(- ))
      end do
      deallocate(cm, Ktt, Knt, D, m, Cu)
      return
      end subroutine

!  for (i in seq(length(x.new)))
!  {
!    cov.new.train <- t(K[train, test[i]])
!    cov.new.new <- K[test[i], test[i]]
!    m <- cov.new.train %*% c.mean
!    var <- cov.new.new - cov.new.train %*% solve(cov.train.train + Dinv) %*% t(cov.new.train)
!    me <- .logistic.normal.integral(m, var)
!    pre <- .sigmoid(.sample.gaussian(m, var))
!    pred.samples <- c(pred.samples , pre)
!    pred.means   <- c(pred.means, me)
!  }
!  m.new <- K[test, train] %*% c.mean
!  cov.new.new <- K[test, test] - K[test, train] %*%
!    solve(cov.train.train + Dinv) %*% K[train, test]
