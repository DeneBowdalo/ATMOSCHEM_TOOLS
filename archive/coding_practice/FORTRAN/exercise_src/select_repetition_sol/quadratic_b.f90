PROGRAM Quadratic_b
   IMPLICIT NONE 
   REAL :: a,b,c, factor, root1, root2
!
! Input a,b and c
!
   WRITE(*,*)'Input a, b, c'
   READ(*,*)a,b,c
!
! Evaluate factor to be fed to sqrt function
!
   factor = b*b - 4.0 * a * c
!
! check that factor is .ge.0 If -ve stop
! with apropriate error message. If factor
! equlas zero we can use the repeated roots message
! else we can evaluate roots as previously
!
   IF(factor.LT.0.0)THEN
       WRITE(*,*)'b*b - 4*a*c is negative'
       STOP
   ELSE IF(factor.EQ.0.0)THEN
       WRITE(*,*)'Repeated roots'
       root1 = (-1.0*b)/(2.0*a)
       root2 = root1
   ELSE 
      factor = SQRT(factor)
      root1 = (-1.0*b + factor)/(2.0*a)
      root2 = (-1.0*b - factor)/(2.0*a)
   END IF
!
! and output results
!
   WRITE(*,*)'The equation has roots ',root1,root2
END PROGRAM Quadratic_b
