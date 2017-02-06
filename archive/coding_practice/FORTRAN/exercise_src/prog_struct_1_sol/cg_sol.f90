Program congugate_gradient

  ! Solve a system of linear equations Ax=b for x where
  ! A is a real symmetric matrix

  Implicit None

  Real, Dimension( :, : ), Allocatable :: a
  Real, Dimension(    : ), Allocatable :: b
  Real, Dimension(    : ), Allocatable :: x

  Real :: tolerance

  Integer :: n
  Integer :: iterations

  Call read_input( n, tolerance )
  Call allocate_arrays( n )
  Call initialize_arrays( a, b )
  Call generate_symm_pos_def( a )
  Call cg_solve( a, b, tolerance, iterations, x )
  Call write_results( a, b, iterations, x )
  Call deallocate_arrays

Contains

  Subroutine read_input( n, tolerance )

    ! Read the input data

    Integer, Intent( Out ) :: n
    Real   , Intent( Out ) :: tolerance

    Real, Parameter :: min_tolerance = 1e-6
    
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

  Subroutine initialize_arrays( a, b )

    ! Set up the arrays with random numbers

    Real, Dimension( :, : ), Intent( Out ) :: a
    Real, Dimension(    : ), Intent( Out ) :: b

    Integer :: n
    Integer :: i

    n = Size( a, Dim = 1 )

    Call Random_number( a )
    Call Random_number( b )

    ! Make sure A is diagonal dominant

    a = a - 0.5

    Do i = 1, n
       a( i, i ) = a( i, i ) + 10
    End Do

    a = a / 10.0

  End Subroutine initialize_arrays

  Subroutine generate_symm_pos_def( a )

    ! From a general matrix generate a real symmetric
    ! positive definite matrix

    Real, Dimension( :, : ), Intent( InOut ) :: a

    a = Matmul( Transpose( a ), a )

  End Subroutine generate_symm_pos_def

  Subroutine cg_solve( a, b, tolerance, iterations, x )

    ! Solve ( to the given tolerance ) the set of equations
    ! Ax=b using the conjugate gradient method. A must be
    ! real, symmetric positive definite. Failure to converge
    ! is indicated by returning the negative of the number
    ! of iterations taken

    Real   , Dimension( :, : ), Intent( In    ) :: a
    Real   , Dimension(    : ), Intent( In    ) :: b
    Real                      , Intent( In    ) :: tolerance
    Integer                   , Intent(   Out ) :: iterations
    Real   , Dimension(    : ), Intent(   Out ) :: x

    Real, Dimension(    : ), Allocatable :: r
    Real, Dimension(    : ), Allocatable :: p

    Real :: alpha
    Real :: beta
    Real :: length_b
    Real :: length_r_squared

    Integer :: n
    Integer :: max_iterations
    Integer :: error

    n = Size( a, Dim = 1 )

    max_iterations = n * n

    Allocate( r( 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of r failed'
       Stop
    End If

    Allocate( p( 1:n ), Stat = error )
    If( error /= 0 ) Then
       Write( *, * ) 'Allocation of p failed'
       Stop
    End If

    length_b = Sqrt( Dot_product( b, b ) )

    x = b

    r = b - Matmul( a, x )
    p = r

    iterations = 0
    Do

       iterations = iterations + 1
       If(  iterations > max_iterations ) Then
          iterations = - iterations
          Exit
       End If

       length_r_squared = Dot_product( r, r )
       alpha = length_r_squared / Dot_product( p, Matmul( a, p ) )

       x = x + alpha * p
       r = r - alpha * Matmul( a, p )

       If( Sqrt( Dot_product( r, r ) ) / length_b <= tolerance ) Then
          Exit
       End If
       
       beta = Dot_product( r, r ) / length_r_squared 
       p = r + beta * p

    End Do

    Deallocate( p )
    Deallocate( r )

  End Subroutine cg_solve

  Subroutine write_results( a, b, iterations, x )

    ! Tell the world what happened

    Real   , Dimension( :, : ), Intent( In ) :: a
    Real   , Dimension(    : ), Intent( In ) :: b
    Integer                   , Intent( In ) :: iterations
    Real   , Dimension(    : ), Intent( In ) :: x

    Real, Dimension(    : ), Allocatable :: r

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

  Subroutine deallocate_arrays

    ! Deallocate the arrays
 
    Deallocate( x )
    Deallocate( b )
    Deallocate( a )

  End Subroutine deallocate_arrays

End Program congugate_gradient
