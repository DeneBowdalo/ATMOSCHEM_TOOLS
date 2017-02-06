PROGRAM Armstrong
  IMPLICIT NONE 
  INTEGER :: a,b,c,integer_number
!
! You may have set the first loop to 1. This is
! a question of is 001 a three digit number?
!
  DO a = 0, 9
     DO b = 0,9
        DO c = 0,9
!
! Note that we need to deal with individual digits of a number
! so we will build up the number from three digits
!
           integer_number = 100*a + 10*b + c
!
! Check for the Armstrong condition
!
           IF(integer_number.EQ.a*a*a+b*b*b+c*c*c)THEN
              WRITE(*,*)integer_number,' is an armstrong number'
           END IF
        END DO
     END DO
  END DO
END PROGRAM Armstrong
