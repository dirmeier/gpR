!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Fortran routines for calculation of binomial class labels
!!!
!!! @author: Simon Dirmeier
!!! @email: simon.dirmeier@gmx.de
!!!

!! compute the binomial class labels
      subroutine predict_binomial(na, ntrain, nnew, K, ctrain, sigtrain, xnew,& 
                                 trainidx, newidx, mnew, ps, pm, Kn)
! load modules
      use gpr_approx
      use gpr_sampler
      use gpr_linalg
      implicit none
      integer, intent(in)           :: na, ntrain, nnew
      integer, intent(in)           :: trainidx(ntrain), newidx(nnew)
      double precision, intent(in)  :: ctrain(ntrain), K(na, na) 
      double precision, intent(in)  :: sigtrain(ntrain), xnew(nnew)
      double precision, intent(out) :: mnew(nnew), ps(nnew), pm(nnew), Kn(nnew, nnew)
! worker arrays
      double precision :: knn, me, pre, var, m, rn
      double precision, allocatable :: cm(:), Ktt(:, :), D(:, :), Knt(:), Cu(:, :)
      integer :: i, cidx, nnewidxst
!
      allocate(cm(ntrain), Ktt(ntrain, ntrain), D(ntrain, ntrain))
      allocate(Knt(ntrain), Cu(ntrain, ntrain))
!
      Ktt = K(1:ntrain, 1:ntrain)
! compute log probability of the class mapping of c.train given
! the posterior mean
! (see references: Barber p406, Equation~19.5.24;
!                  Rasmussen p44, Equation~3.15 & Equation~3.21)
      cm = ctrain - sigtrain
      D = 0.0
! worker       
      Cu  = Ktt + D
      call solve(ntrain, Cu)
! compute the inverse of the Hessian of the log-likelihood      
      forall(i = 1:ntrain) D(i, i) = sigtrain(i) * (1 - sigtrain(i))
      do i = 1, nnew
          cidx = newidx(i)
! get the covariance of the training points and a new input          
          Knt = K(1:ntrain, cidx)
! get the variance of the new input point          
          knn = K(cidx, cidx)
! compute the posterior mean of the predictive distribution
! (see references: Barber p406, Equation~19.5.24;
!                  Rasmussen p44, Equation~3.21)          
          m   = DOT_PRODUCT(Knt, cm)
! compute the posterior variance of the predictive distribution
! (see references: Barber p406, Equation~19.5.26;
!                  Rasmussen p44, Equation~3.24)          
          var = knn - DOT_PRODUCT(MATMUL(Knt, Cu), Knt)
! compute the mean class probability mapping
! (see references: Barber p406, Equation~19.5.27;
!                  Rasmussen p44, Equation~3.25)          
          me  = lognint(m, var)
          rn = rnorm(m, var)
! compute a sample from the posterior class mapping          
          pre =  1 / (1 + exp(-1 * rn ))
          ps(i) = pre
          pm(i) = me
      end do
      nnewidxst = ntrain + 1
! mean of posterior       
      mnew = MATMUL(K(nnewidxst:na, 1:ntrain), cm)
! variance of posterior      
      Kn = K(newidx, newidx) / MATMUL(MATMUL(K(nnewidxst:na, 1:ntrain), Cu), &
                                      K(1:ntrain, nnewidxst:na))
!
      deallocate(cm, Ktt, Knt, D, Cu)
      return
      end subroutine

