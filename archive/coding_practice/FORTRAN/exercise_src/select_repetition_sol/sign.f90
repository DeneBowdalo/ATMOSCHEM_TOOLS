PROGRAM posneg
  IMPLICIT NONE
  INTEGER :: input

  DO !infinite do loop
     WRITE (*,*) 'ENTER an integer (or 99999 to finish)'
     READ (*,*) input
!
!test for stopping criteria
!
     IF (input == 99999)THEN
        EXIT  !exit the do loop
     END IF
!
!test the input
!
     IF (input > 0) THEN
      WRITE(*,*) 'Positive'
     ELSE IF (input == 0) THEN
      WRITE(*,*) 'Zero'
     ELSE
       WRITE(*,*)'Negative'
     END IF
  END DO
  
END PROGRAM posneg
