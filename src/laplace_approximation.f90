      subroutine laplace_approximation(n, c, K, y, sig, D)
      implicit none
      integer, intent(in) :: n
      double precision, intent(in) :: c(n), K(n, n)
      double precision, intent(out) :: y(n), sig(n), D(n, n)

      integer :: i
      do i = 1, n
          y(i) = 1.0
      end do
      
      sig(:) = 2.0
      D(:, :) = 3.0
      print *, n
      print *, sig
      return
      end subroutine

