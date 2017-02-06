MODULE Sort_it
  ! A module to sort arrays of reals or integers
  Use numbers_module
  IMPLICIT NONE 
  PRIVATE
  PUBLIC sort
!
!an interface to overload a sorting rountine
!
  INTERFACE sort
     MODULE PROCEDURE sort_int
     MODULE PROCEDURE sort_real
  END INTERFACE

CONTAINS

  SUBROUTINE sort_int( a )

    ! Integer sorting routine. Method is selection sort.

    INTEGER, DIMENSION( : ), INTENT( InOut ) :: a

    INTEGER, DIMENSION(1) :: location

    INTEGER :: n
    INTEGER :: swap
    INTEGER :: i, j

    n = SIZE( a )

    DO i = 1, n - 1
!
! Note two things; firstly minloc returns an array.
! Secondly in returns the location of the minimum within the
! array segment passed to it, not within the whole array. Hence
! we need to add (i-1) to its value to take account of this offset
!
       location = MINLOC( a( i:n ) )
       j = location( 1 ) + i - 1
       swap   = a( j )
       a( j ) = a( i )
       a( i ) = swap

    END DO

  END SUBROUTINE sort_int

  SUBROUTINE sort_real( a )

    ! Real sorting routine. Method is selection sort.
    
    REAL( wp ), DIMENSION( : ), INTENT( InOut ) :: a
    
    INTEGER, DIMENSION(1) :: location
    
    INTEGER :: n
    REAL( wp ) :: swap
    INTEGER :: i, j
    
    n = SIZE( a )
    
    DO i = 1, n - 1
!
! Note two things; firstly minloc returns an array.
! Secondly in returns the location of the minimum within the
! array segment passed to it, not within the whole array. Hence
! we need to add (i-1) to its value to take account of this offset
!
       location = MINLOC( a( i:n ) )
       j = location( 1 ) + i - 1
       swap   = a( j )
       a( j ) = a( i )
       a( i ) = swap
       
    END DO
    
  END SUBROUTINE sort_real
  
END MODULE Sort_it
    
