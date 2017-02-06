PROGRAM Calc_pi
   IMPLICIT NONE 
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,70) 
   INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6,35)
   INTEGER, PARAMETER :: prec = dp !change to sp for single precision
   REAL(KIND = prec) :: a,b,p, new_p 
   INTEGER :: i=0
   
!
!Set the initial values
!
   a = sqrt(2.0_prec)
   b = 0.0_prec
   p = 2.0_prec+a
!
!iterate until two subsequent values of p are equal
!
   do 
      b = sqrt(a)*(1.0_prec+b)/(a+b)
      a = (sqrt(a) + 1.0_prec /sqrt(a))/2.0_prec
      new_p = p*b*(1.0_prec+a)/(1.0_prec+b)
      i = i + 1
      if(new_p-p==0.0_prec)exit
      p = new_p
   end do
!
! write out our value of pi
!
   write(*,*) p
   write(*,*) 'number of iterations ', i
END PROGRAM Calc_pi
