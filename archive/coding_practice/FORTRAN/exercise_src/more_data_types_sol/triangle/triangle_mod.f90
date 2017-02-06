MODULE TRIANGLE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: triang, Cos_or_pyth
!
! derived type for a triangle
!
  TYPE triang
     REAL :: A
     REAL :: B
     REAL :: theta
     LOGICAL :: rightang
  END TYPE triang


CONTAINS


  SUBROUTINE Cos_or_pyth(tri,C)
    IMPLICIT NONE
    TYPE(triang), INTENT(IN) :: tri 
    REAL, INTENT(OUT) :: C
    INTEGER, PARAMETER :: cosine = 1
    INTEGER, PARAMETER :: pythag = 2
    INTEGER :: method_to_use
    
!
!check if there is an angle included then 
!pick the appropriate method to calculate the length of the 3rd side
!
    IF(tri%rightang) THEN
       method_to_use = pythag
    ELSE
       method_to_use = cosine
    END IF
    
    SELECT CASE (method_to_use)
    CASE(cosine)  
       C = SQRT(tri%A**2.0+tri%B**2.0-2.0*tri%A*tri%B*cos(tri%theta)) !cosine rule
    CASE(pythag)
       C = SQRT(tri%A**2.0+tri%B**2.0) !Pythagoras' Theorem
    CASE DEFAULT
       WRITE(*,*) 'Unknown method'
       STOP
    END SELECT
    
    
  END SUBROUTINE Cos_or_pyth
  
  
END MODULE TRIANGLE
