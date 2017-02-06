PROGRAM Cosine
  IMPLICIT NONE
  REAL :: A, B, C, theta

!
!read in the sides and the angle  
! 
  WRITE(*,*) 'How long is the first side?'
  READ(*,*) A
  WRITE(*,*) 'How long is the second side?'
  READ(*,*) B
  WRITE(*,*)'What is the angle between them (in radians)?'
  READ(*,*) theta
!
!use the cosine rule to calculate the length of the 3rd side
! 
  C = SQRT(A**2.0+B**2.0-2.0*A*B*cos(theta))
  WRITE(*,*) 'The third side is', C


END PROGRAM Cosine
