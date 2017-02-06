Module base_data

  Use numbers_module

  Implicit None

  Real( wp ), Dimension( :, : ), Allocatable, Public :: a
  Real( wp ), Dimension(    : ), Allocatable, Public :: b
  Real( wp ), Dimension(    : ), Allocatable, Public :: x

  Public :: allocate_arrays
  Public :: deallocate_arrays

  Private

Contains

  Subroutine allocate_arrays( n )

    ! Allocate the data

    Integer, Intent( In ) :: n

    Integer :: error

    Allocate( a( 1:n, 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of a failed'
       Stop
    End If

    Allocate( b( 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of b failed'
       Stop
    End If

    Allocate( x( 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of x failed'
       Stop
    End If

  End Subroutine allocate_arrays

  Subroutine deallocate_arrays

    ! Deallocate the arrays
 
    Deallocate( x )
    Deallocate( b )
    Deallocate( a )

  End Subroutine deallocate_arrays

End Module base_data
