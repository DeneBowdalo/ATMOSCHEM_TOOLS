PROGRAM TRI_MAIN
  USE TRIANGLE
  USE CONVERT, ONLY : deg_to_rad 
  IMPLICIT NONE
  REAL :: C
  TYPE(triang) :: tri
  LOGICAL :: deg

  WRITE(*,*)'Is it a right angled triangle? (T or F) '
  READ(*,*) tri%rightang

! 
!if the locical is true then it is a right angled triangle
!so we don't need to ask about the angle nor include it in
!the subroutine call
!
  IF(tri%rightang) THEN    
     WRITE(*,*) 'How long is the first side?'
     READ(*,*) tri%A
     WRITE(*,*) 'How long is the second side?'
     READ(*,*) tri%B
     tri%theta  = 2.0*ATAN(1.0)
     CALL Cos_or_pyth(tri,C)
  ELSE
     WRITE(*,*) 'How long is the first side?'
     READ(*,*) tri%A
     WRITE(*,*) 'How long is the second side?'
     READ(*,*) tri%B
     WRITE(*,*)'What is the angle between them?'
     READ(*,*) tri%theta
     WRITE(*,*) 'Is that in Degrees? (T or F)'
     READ(*,*) deg
     IF(deg) CALL deg_to_rad(tri%theta)
     CALL Cos_or_pyth(tri,C)
  END IF
  WRITE(*,*)'The third side is ', C

END PROGRAM TRI_MAIN
