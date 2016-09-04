      subroutine descent(b, X, Y, n, p)
      integer, intent(in) :: n, p 
      real, intent(out) :: b(p)
      real, intent(in) :: X(n, p), Y(n)
      
      real :: Xtran(p, n), TXX(p, p), TXY(p), bold(p), lamb, softthresh, maxim, si
      integer :: i, j, EXL(p), betaIdxMin(p-1)
      lamb = 0.0       
      EXL = (/ (i,i=1,p) /)      
      Xtran = TRANSPOSE(X)
      TXX = MATMUL(Xtran, X)      
      TXY = MATMUL(Xtran, Y)      
      b = 0.0
      bold = 1.0      
      do while (sum(abs(b-bold)) > 0.001)      
       print *, b, " -> ",  bold 
       bold = b
       do i = 1, p
        betaIdxMin = pack(EXL, EXL /= i)        
        softthresh = TXY(i)
        do j = 1, (p-1)
         softthresh = softthresh - TXX(i, betaIdxMin(j)) * b(betaIdxMin(j))            
        end do
        si = merge(1.0, merge(0.0, -1.0, softthresh == 0.0), softthresh > 0.0)         
        maxim = merge(abs(softthresh - lamb), 0.0, abs(softthresh - lamb) > 0)
        b(i) = (si * maxim) / TXX(i,i)             
       end do      
      end do     
      end subroutine 
