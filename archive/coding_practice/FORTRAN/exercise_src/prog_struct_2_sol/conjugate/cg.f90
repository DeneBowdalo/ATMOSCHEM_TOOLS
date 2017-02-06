Module cg

  Use numbers_module

  Implicit None

  Public :: cg_solve

  Private

Contains

  Subroutine cg_solve( a, b, tolerance, iterations, x )

    ! Solve ( to the given tolerance ) the set of equations
    ! Ax=b using the conjugate gradient method. A must be
    ! real, symmetric positive definite. Failure to converge
    ! is indicated by returning the negative of the number
    ! of iterations taken

    Real( wp )   , Dimension( :, : ), Intent( In    ) :: a
    Real( wp )   , Dimension(    : ), Intent( In    ) :: b
    Real( wp )                      , Intent( In    ) :: tolerance
    Integer                         , Intent(   Out ) :: iterations
    Real( wp )   , Dimension(    : ), Intent(   Out ) :: x

    Real( wp ), Dimension(    : ), Allocatable :: r
    Real( wp ), Dimension(    : ), Allocatable :: p

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

End Module cg
