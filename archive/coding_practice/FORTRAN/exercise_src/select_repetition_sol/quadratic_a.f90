PROGRAM Quadratic_a
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
! with apropriate error message
!
   IF(factor.LT.0.0)THEN
       WRITE(*,*)'b*b - 4*a*c is negative'
       STOP
   END IF
!
! As factor is ge.0 we can take its sqrt
! and evaluate roots
!
   factor = SQRT(factor)
   root1 = (-1.0*b + factor)/(2.0*a)
   root2 = (-1.0*b - factor)/(2.0*a)
!
! and output results
!
   WRITE(*,*)'The equation has roots ',root1,root2
END PROGRAM Quadratic_a
