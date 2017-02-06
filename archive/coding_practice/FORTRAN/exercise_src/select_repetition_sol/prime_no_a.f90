PROGRAM Prime_a
  IMPLICIT NONE 
  INTEGER :: input_int,i
  LOGICAL :: prime
!
! Set a flag to indicate if we have found a prime
! Default to .true.
!
  prime = .TRUE.
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
! test - try and work out the logic of the next line
!
  DO i = 2, INT(SQRT(REAL(input_int)))+1
     IF(MOD(input_int,i).EQ.0)THEN
        WRITE(*,*)'Number ',input_int,' is not prime'
        WRITE(*,*)'One root is',i
        prime = .FALSE.
        EXIT
     END IF
  END DO
!
! THE preceding loop will have set prime to .false. if any number
! between 2 and the sqrt(i) is a factor. So if prime is still true
! then the number is prime so output a suitable message
!
  IF(prime)THEN
     WRITE(*,*)'Number ',input_int,'is prime'
  END IF
END PROGRAM Prime_a
