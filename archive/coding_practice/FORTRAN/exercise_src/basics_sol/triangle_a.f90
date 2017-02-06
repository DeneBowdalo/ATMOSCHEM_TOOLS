PROGRAM Pythag
  IMPLICIT NONE
  REAL :: A, B, C
!
!read in the length of the sides
!  
  WRITE(*,*) 'How long is the 1st side?'
  READ(*,*) A
  WRITE(*,*) 'How long is the 2nd side?'
  READ(*,*) B
!
!use Pythagoras' theorem to calculate the length of the 3rd side
!
  C = SQRT(A**2.0+B**2.0)
  WRITE(*,*) 'The hypotenuse is', C


END PROGRAM Pythag
