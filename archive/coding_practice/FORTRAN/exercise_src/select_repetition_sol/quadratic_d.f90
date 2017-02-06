PROGRAM Quadratic_d
  IMPLICIT NONE 
  REAL    :: a,b,c, factor
  COMPLEX :: root1, root2
  COMPLEX, PARAMETER :: runity = (1.0,0.0)
  !
  ! Note we have defined a new constant
  ! called r(eal)unity. Input a,b and c
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
     ! Construct factor and check for -ve,0,+ve. Note that if root1
     ! is of the form (a,0) (e.g. real) then it will be written as
     ! (a,0). We can cast its type to real so it prints as a
     !
     factor = b*b - 4.0 * a * c
     IF(factor.EQ.0.0)THEN
        root1 = (-runity*b)/(2.0*a*runity)
        WRITE(*,*)'Repeated root; the solution is', REAL(root1), ' twice'
     ELSE 
        !
        ! We are now solving using a numerically stable formula. As root1 and root2 are complex
        ! we do not need to check for factor less than zero before feeding to sqrt
        ! As factor is of type real the usual trick of factor = sqrt(factor) would
        ! be wrong.
        !	
        IF(b.LT.0)THEN
           root1 = (-1.0*b*runity + SQRT(factor*runity))/(2.0*a)
        ELSE
           root1 = (-1.0*b*runity - SQRT(factor*runity))/(2.0*a)
        END IF

        root2 = c/a/root1

        !
        ! However we will check before output so that if the roots are real we
        ! can cast the complex numbers to type real
        !
        IF(factor.LT.0.0)THEN
           WRITE(*,*)'The equation has roots ',root1,root2
        ELSE
           WRITE(*,*)'The equation has roots ',REAL(root1),REAL(root2)
        END IF
     END IF
  END IF
END PROGRAM Quadratic_d
