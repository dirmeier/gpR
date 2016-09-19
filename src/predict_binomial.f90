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
      double precision function lognint(m, var) result(e)
      implicit none
      double precision, intent(in) :: m, var
      double precision, parameter :: PI = 3.1415926535897932
      double precision            :: lam
      lam = sqrt(PI / 4)
      e = .5 + .5 * erf(lam * m / sqrt(1.0 + (2 * (lam ** 2.0) * var)))
      end function lognint

! draw from a Gaussian distribution using polar Box-Muller transform
      double precision function rnorm(m, sd)
      implicit none
      double precision, intent(in) :: m, sd
      double precision :: x, y, tr

      tr = 2.0
      do while(tr >= 1)
          x = 2.0 * rand() - 1.0
          y = 2.0 * rand() - 1.0
          tr = x * x + y * y
      end do
      tr = sqrt((-2.0 * log(tr)) / tr)
      rnorm = tr * x
      end function rnorm


! compute the binomial class labels
      subroutine predict_binomial(na, ntrain, nnew, K, ctrain, sigtrain, xnew, trainidx, newidx, mnew, ps, pm, Kn)
      implicit none
      integer, intent(in)           :: na, ntrain, nnew
      integer, intent(in)           :: trainidx(ntrain), newidx(nnew)
      double precision, intent(in)  :: ctrain(ntrain), K(na, na) 
      double precision, intent(in)   :: sigtrain(ntrain), xnew(nnew)
      double precision, intent(out) :: mnew(nnew), ps(nnew), pm(nnew), Kn(nnew, nnew)
!
      double precision :: knn, me, pre, var, m, rn
      double precision, allocatable :: cm(:), Ktt(:, :), D(:, :), Knt(:, :), Cu(:, :)
      integer :: i, cidx
!
      allocate(cm(ntrain), Ktt(ntrain, ntrain), D(ntrain, ntrain))
      allocate(Knt(ntrain, 1), Cu(ntrain, ntrain))
!
      Ktt = K(1:ntrain, 1:ntrain)
      cm = ctrain - sigtrain
      D = 0.0
      Cu  = Ktt + D
      call solve(ntrain, Cu)
      forall(i = 1:ntrain) D(i, i) = sigtrain(i) * (1 - sigtrain(i))
      do i = 1, nnew
          cidx = newidx(i)
          Knt = K(1:trainidx, newidx(i))
          knn = K(cidx, cidx)
          m   = sum(MATMUL(Knt, cm))
          var = knn - sum(MATMUL(MATMUL(TRANSPOSE(Knt), Cu), Knt))
          me  = lognint(m, var)
          rn = rnorm(m, var)
          pre =  1 / (1 + exp(-1 * rn ))
          ps(i) = pre
          pm(i) = me
      end do
      mnew = MATMUL( K((ntrain+1):nnew, 1:ntrain), cm)
      Kn = K(newidx, newidx) / MATMUL(MATMUL(K((ntrain+1):nnew, 1:ntrain), Cu), K((ntrain+1):nnew, 1:ntrain))

      deallocate(cm, Ktt, Knt, D, Cu)
      return
      end subroutine

