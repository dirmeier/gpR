!!
!!
!
      module gpr_util
      implicit none
!
      contains 
!
! calculate the logistic normal integral using the error function
      double precision function lognint(m, var) result(e)
      implicit none
      double precision, intent(in) :: m, var
      double precision, parameter :: PI = 3.1415926535897932
      double precision            :: lam
      lam = sqrt(PI / 4)
      e = DBLE(.5 + .5 * erf(lam * m / sqrt(1.0 + (2 * (lam ** 2.0) * var))))
      end function lognint
!
!
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
!
!
! draw from a Gaussian distribution using polar Box-Muller transform
      double precision function rnorm(m, sd)
      implicit none
      double precision, intent(in) :: m, sd
      double precision :: x, y, tr
      real :: r
      tr = 2.0
      do while(tr >= 1)
          call random_number(r)
          x = 2.0 * r - 1.0
          call random_number(r)
          y = 2.0 * r - 1.0
          tr = x * x + y * y
      end do
      tr = sqrt((-2.0 * log(tr)) / tr)
      rnorm = tr * x
      end function rnorm
!      
!
      end module gpr_util

