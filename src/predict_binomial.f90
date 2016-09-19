!! Fortran routines for calculation of binomial class labels
!!
! compute the binomial class labels
      subroutine predict_binomial(na, ntrain, nnew, K, ctrain, sigtrain, xnew, trainidx, newidx, mnew, ps, pm, Kn)
      use gpr_util
      implicit none
      integer, intent(in)           :: na, ntrain, nnew
      integer, intent(in)           :: trainidx(ntrain), newidx(nnew)
      double precision, intent(in)  :: ctrain(ntrain), K(na, na) 
      double precision, intent(in)   :: sigtrain(ntrain), xnew(nnew)
      double precision, intent(out) :: mnew(nnew), ps(nnew), pm(nnew), Kn(nnew, nnew)
!
      double precision :: knn, me, pre, var, m, rn
      double precision, allocatable :: cm(:), Ktt(:, :), D(:, :), Knt(:), Cu(:, :)
      integer :: i, cidx
!
      allocate(cm(ntrain), Ktt(ntrain, ntrain), D(ntrain, ntrain))
      allocate(Knt(ntrain), Cu(ntrain, ntrain))
!
      Ktt = K(1:ntrain, 1:ntrain)
      cm = ctrain - sigtrain
      D = 0.0
      Cu  = Ktt + D
      call solve(ntrain, Cu)
      forall(i = 1:ntrain) D(i, i) = sigtrain(i) * (1 - sigtrain(i))
      do i = 1, nnew
          cidx = newidx(i)
          Knt = K(1:ntrain, newidx(i))
          knn = K(cidx, cidx)
          m   = DOT_PRODUCT(Knt, cm)
          var = knn - DOT_PRODUCT(MATMUL(Knt, Cu), Knt)
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

