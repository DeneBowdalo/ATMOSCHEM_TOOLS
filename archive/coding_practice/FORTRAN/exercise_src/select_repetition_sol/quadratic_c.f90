PROGRAM Quadratic_c
   IMPLICIT NONE 
   REAL :: a,b,c, factor, root1, root2
!
! Input a,b and c
!
   WRITE(*,*)'Input a, b, c'
   READ(*,*)a,b,c
!
! First check status of a
!
   IF(a.EQ.0.0)THEN
!
! If a is zero then check b
!
       IF(b.NE.0.0)THEN
           WRITE(*,*)'Linear equation: Root is',-1.0*c/b
       ELSE
!
! and check c
!
           IF(c.EQ.0.0)THEN
               WRITE(*,*)'Solution undefined. All parameters = 0.0'
           ELSE
               WRITE(*,*)'No solution'
           END IF
       END IF
   ELSE
!
! We have a quadratic equation so we are back on familiar ground
! Construct factor and check for -ve,0,+ve
!
      factor = b*b - 4.0 * a * c
      IF(factor.LT.0.0)THEN
         WRITE(*,*)'b*b - 4*a*c is negative'
         STOP
      ELSE IF(factor.EQ.0.0)THEN
         root1 = (-1.0*b)/(2.0*a)
         WRITE(*,*)'Repaeted root; the solution is', root1, ' twice'
      ELSE 
         factor = SQRT(factor)
         root1 = (-1.0*b + factor)/(2.0*a)
         root2 = (-1.0*b - factor)/(2.0*a)
          WRITE(*,*)'The equation has roots ',root1,root2
      END IF
   END IF
END PROGRAM Quadratic_c
