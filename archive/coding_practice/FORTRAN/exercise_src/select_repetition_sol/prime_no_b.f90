PROGRAM Prime_b
  IMPLICIT NONE 
  INTEGER :: input_int,i,j
  LOGICAL :: prime
!
! Now we need a general do loop for reading. This is
! a simple error check for the input number being > 2
!
  DO
     WRITE(*,*)'Input an integer number > 2'
     READ(*,*)input_int
     IF(input_int.LE.2)THEN
        WRITE(*,*)'Error. Number must be > 2'
     ELSE
        EXIT
     END IF
  END DO
!  
  WRITE(*,*)1
  WRITE(*,*)2
  DO j = 3, input_int
     prime = .TRUE.
     DO i = 2, INT(SQRT(REAL(j)))+1
        IF(MOD(j,i).EQ.0)THEN
           prime = .FALSE.
           EXIT
        END IF
     END DO
     IF(prime)THEN
        WRITE(*,*)j
     END IF
  END DO
END PROGRAM Prime_b
