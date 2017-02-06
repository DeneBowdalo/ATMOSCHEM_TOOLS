PROGRAM TRI_MAIN
  Use numbers_module
  USE TRIANGLE
  USE CONVERT, ONLY : deg_to_rad 
  IMPLICIT NONE
  REAL( wp ) :: A, B, theta, C
  LOGICAL :: angle, deg

  WRITE(*,*)'Is it a right angled triangle? (T or F) '
  READ(*,*) angle

! 
!if the locical is true then it is a right angled triangle
!so we don't need to ask about the angle nor include it in
!the subroutine call
!
  IF(angle) THEN    
     WRITE(*,*) 'How long is the first side?'
     READ(*,*) A
     WRITE(*,*) 'How long is the second side?'
     READ(*,*) B
     CALL Cos_or_pyth(A,B,C)
  ELSE
     WRITE(*,*) 'How long is the first side?'
     READ(*,*) A
     WRITE(*,*) 'How long is the second side?'
     READ(*,*) B
     WRITE(*,*)'What is the angle between them?'
     READ(*,*) theta
     WRITE(*,*) 'Is that in Degrees? (T or F)'
     READ(*,*) deg
     IF(deg) CALL deg_to_rad(theta)
     CALL Cos_or_pyth(A,B,C,theta)
  END IF
  WRITE(*,*)'The third side is ', C

END PROGRAM TRI_MAIN
