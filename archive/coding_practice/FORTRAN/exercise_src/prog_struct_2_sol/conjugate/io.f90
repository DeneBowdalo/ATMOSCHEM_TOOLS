Module io

  Use numbers_module

  Implicit None

  Public :: read_input
  Public :: write_results

  Private

Contains

  Subroutine read_input( n, tolerance )

    ! Read the input data

    Integer   , Intent( Out ) :: n
    Real( wp ), Intent( Out ) :: tolerance

    Real( wp ), Parameter :: min_tolerance = 1e-6_wp
    
    Write( *, * ) 'How big should the problem be ?'
    Read ( *, * ) n

    Write( *, * ) 'And the tolerance ?'
    Do
       Read ( *, * ) tolerance
       If( tolerance >= min_tolerance ) Then
          Exit
       Else
          Write( *, * ) 'Tolerance smaller than the minimum allowed: ', &
               min_tolerance
          Write( *, * ) 'Please input a larger value'
       End If
    End Do

  End Subroutine read_input

  Subroutine write_results( a, b, iterations, x )

    ! Tell the world what happened

    Real( wp )   , Dimension( :, : ), Intent( In ) :: a
    Real( wp )   , Dimension(    : ), Intent( In ) :: b
    Integer                         , Intent( In ) :: iterations
    Real( wp )   , Dimension(    : ), Intent( In ) :: x

    Real( wp ), Dimension(    : ), Allocatable :: r

    Integer :: n
    Integer :: error
    Integer :: i
      
    n = Size( a, Dim = 1 )

    Allocate( r( 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of r failed'
       Stop
    End If

    Write( *, * ) 'Input Matrix'
    Write( *, * ) '------------'
    Do i = 1, n
       Write( *, '( 1x, 100( f5.2, 1x ) )' ) a( i, : )
    End Do

    Write( *, * ) 
    Write( *, * ) 'Input RHS Vector'
    Write( *, * ) '----------------'
    Do i = 1, n
       Write( *, '( 1x, f5.2 )' ) b( i )
    End Do
    Write( *, * )

    If( iterations > 0 ) Then

       Write( *, *  ) 'The CG method converged in ', &
            Abs( iterations ), ' iterations.'

       Write( *, * ) 
       Write( *, * ) 'Solution'
       Write( *, * ) '--------'
       Do i = 1, n
          Write( *, '( 1x, f8.5 )' ) x( i )
       End Do
       
       r = b - Matmul( a, x )
       
       Write( *, * ) 
       Write( *, * ) 'Error vector'
       Write( *, * ) '------------'
       Do i = 1, n
          Write( *, '( 1x, f8.5 )' ) r( i )
       End Do

    Else

       Write( *, *  ) 'The CG method failed to converge in ', &
            Abs( iterations ), ' iterations.'

    End If

    Deallocate( r )

  End Subroutine write_results

End Module io
