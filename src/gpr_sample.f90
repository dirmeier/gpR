!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Module for sampling from distributions.
!!!
!!! @author: Simon Dirmeier
!!! @email: simon.dirmeier@gmx.de
!!!
      module gpr_sampler
      implicit none
!! code block
      contains
!! draw from a Gaussian distribution using iterative polar Box-Muller transform
      double precision function rnorm(m, sd)
      implicit none
      double precision, intent(in) :: m, sd
      double precision :: x, y, tr
      real :: r
      tr = 2.0
! do iterative Box Muller for numerical stability
      do while(tr >= 1)
          call random_number(r)
          x = 2.0 * r - 1.0
          call random_number(r)
          y = 2.0 * r - 1.0
          tr = x * x + y * y
      end do
      tr = sqrt((-2.0 * log(tr)) / tr)
      rnorm = (tr * x) * sd + m
      end function rnorm
! end of module :(
      end module gpr_sampler

