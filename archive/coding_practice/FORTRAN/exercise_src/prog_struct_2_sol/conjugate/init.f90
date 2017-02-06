Module init

  Use numbers_module
   
  Implicit None

  Public :: initialize_arrays
  Public :: generate_symm_pos_def

  Private

Contains

  Subroutine initialize_arrays( a, b )

    ! Set up the arrays with random numbers

    Real( wp ), Dimension( :, : ), Intent( Out ) :: a
    Real( wp ), Dimension(    : ), Intent( Out ) :: b

    Integer :: n
    Integer :: i

    n = Size( a, Dim = 1 )

    Call Random_number( a )
    Call Random_number( b )

    ! Make sure A is diagonal dominant

    a = a - 0.5_wp

    Do i = 1, n
       a( i, i ) = a( i, i ) + n
    End Do

    a = a / n

  End Subroutine initialize_arrays

  Subroutine generate_symm_pos_def( a )

    ! From a general matrix generate a real symmetric
    ! positive definite matrix

    Real( wp ), Dimension( :, : ), Intent( InOut ) :: a

    a = Matmul( Transpose( a ), a )

  End Subroutine generate_symm_pos_def

End Module init
