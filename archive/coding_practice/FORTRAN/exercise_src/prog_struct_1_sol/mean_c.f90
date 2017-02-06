Program Means

  ! Program to read in data and calculate the arithmetic, harmonic
  ! and geometric means of that data.

  IMPLICIT NONE

  REAL, DIMENSION( : ), ALLOCATABLE :: data

  REAL :: arithmetic_mean
  REAL :: geometric_mean
  REAL :: harmonic_mean

  ! Read in the data
  CALL read_data

  ! Calculate the means
  CALL calc_arithmetic_mean( data, arithmetic_mean )
  CALL calc_geometric_mean ( data, geometric_mean  )
  CALL calc_harmonic_mean  ( data, harmonic_mean   )

  ! Write out the results
  CALL write_results( arithmetic_mean, geometric_mean, &
       harmonic_mean )

  ! Free the memory
  CALL deallocate_data

CONTAINS

  SUBROUTINE read_data

    INTEGER :: n_data
    INTEGER :: error
    INTEGER :: i

    WRITE( *, * ) 'How many data points do you have ?'
    READ ( *, * ) n_data

    ALLOCATE( data( 1:n_data ), Stat = error )
    IF( error /= 0 ) THEN
       WRITE( *, * ) 'Failed to allocate the array data'
    END IF

    DO i = 1, n_data
       WRITE( *, * ) 'Please enter data point ', i
       READ ( *, * ) data( i )
    END DO

  END SUBROUTINE read_data

  SUBROUTINE calc_arithmetic_mean( data, mean )

    REAL, DIMENSION( : ), INTENT( IN    ) :: data
    REAL                , INTENT(   OUT ) :: mean

    mean = SUM( data ) / SIZE( data )

  END SUBROUTINE calc_arithmetic_mean

  SUBROUTINE calc_geometric_mean( data, mean )

    REAL, DIMENSION( : ), INTENT( IN    ) :: data
    REAL                , INTENT(   OUT ) :: mean

    mean = PRODUCT( data ) ** ( 1.0 / REAL( SIZE( data ) ) )

  END SUBROUTINE calc_geometric_mean

  SUBROUTINE calc_harmonic_mean( data, mean )

    REAL, DIMENSION( : ), INTENT( IN    ) :: data
    REAL                , INTENT(   OUT ) :: mean

    mean = REAL( SIZE( data ) ) / SUM( 1.0 / data ) 

  END SUBROUTINE calc_harmonic_mean

  SUBROUTINE write_results( arithmetic_mean, geometric_mean, &
       harmonic_mean )

    REAL, INTENT( IN ) :: arithmetic_mean
    REAL, INTENT( IN ) :: geometric_mean
    REAL, INTENT( IN ) :: harmonic_mean

    WRITE( *, '( 1x, "The arithmetic mean is: ", f10.5 )' ) arithmetic_mean
    WRITE( *, '( 1x, "The geometric  mean is: ", f10.5 )' ) geometric_mean
    WRITE( *, '( 1x, "The harmonic   mean is: ", f10.5 )' ) harmonic_mean

  END SUBROUTINE write_results

  SUBROUTINE deallocate_data

    DEALLOCATE( data )

  END SUBROUTINE deallocate_data

END PROGRAM means
