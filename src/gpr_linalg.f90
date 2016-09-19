!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Module for linear algebra routines.
!!!
!!! @author: Simon Dirmeier
!!! @email: simon.dirmeier@gmx.de
!!!
      module gpr_linalg
      implicit none
!! code block
      contains
!! invert a matrix
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

      end module gpr_linalg

