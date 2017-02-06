MODULE TRIANGLE
  Use numbers_module
  IMPLICIT NONE

CONTAINS


  SUBROUTINE Cos_or_pyth(A,B,C,theta)
    IMPLICIT NONE
    REAL( wp ), INTENT(IN) :: A, B 
    REAL( wp ), INTENT(IN), OPTIONAL :: theta
    REAL( wp ), INTENT(OUT) :: C
    INTEGER, PARAMETER :: cosine = 1
    INTEGER, PARAMETER :: pythag = 2
    INTEGER :: method_to_use
    
!
!check if there is an angle included then 
!pick the appropriate method to calculate the length of the 3rd side
!
    IF( PRESENT(theta)) THEN
       method_to_use = cosine
    ELSE
       method_to_use = pythag
    END IF

    SELECT CASE (method_to_use)
    CASE(cosine)  
       C = SQRT(A**2.0_wp+B**2.0_wp-2.0_wp*A*B*cos(theta)) !cosine rule
    CASE(pythag)
       C = SQRT(A**2.0_wp+B**2.0_wp) !Pythagoras' Theorem
    CASE DEFAULT
       WRITE(*,*) 'Unknown method'
       STOP
    END SELECT
    
    
  END SUBROUTINE Cos_or_pyth
  
  
END MODULE TRIANGLE
