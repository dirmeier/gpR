!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Utility module for approximation algorithms.
!!!
!!! @author: Simon Dirmeier
!!! @email: simon.dirmeier@gmx.de
!!!
      module gpr_approx
      implicit none
!! code block
      contains
!!
! calculate the logistic normal integral using the error function
!!
      double precision function lognint(m, var) result(e)
      implicit none
      double precision, intent(in) :: m, var
      double precision, parameter :: PI = 3.1415926535897932
      double precision            :: lam
      lam = sqrt(PI) / 4
      ! approximation of logistic normal integral using the Gaussian error function
      e = DBLE(.5 + .5 * erf(lam * m / sqrt(1.0 + (2 * (lam ** 2.0) * var))))
      end function lognint
! end of module :(
      end module gpr_approx

