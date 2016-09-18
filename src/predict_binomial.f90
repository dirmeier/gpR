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
      pure function rnorm(m, sd)
      implicit none
      double precision, intent(in) :: m, sd
      double precision             :: x, y, t, rnorm
!
      t = 2.0
      do while(t >= 1)
          x = 2.0 * rand() - 1.0
          y = 2.0 * rand() - 1.0
          t = x * x + y * y
      end do
      t = sqrt((-2.0 * log(t)) / t)
      rnorm = t * x
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
          pre =  1 / (1 + exp(-rnorm(m, var)))
          ps(i) = pre
          pm(i) = me
      end do
i     
      deallocate(cm, Ktt, Knt, D, m, Cu)
      return
      end subroutine

