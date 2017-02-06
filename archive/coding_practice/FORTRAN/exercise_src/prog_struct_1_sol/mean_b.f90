Program means

  ! Program to calculate the arithmetic, geometric
  ! and harmonic average of an arbritary number of
  ! numbers provided by the user.

  Implicit None

  Real, Dimension( : ), Allocatable :: data

  Real :: arithmetic_mean
  Real :: geometric_mean
  Real :: harmonic_mean

  Integer :: n_data
  Integer :: error
  Integer :: i

  Write( *, * ) 'How many data points do you have ?'
  Read ( *, * ) n_data

  Allocate( data( 1:n_data ), Stat = error )
  If( error /= 0 ) Then
     Write( *, * ) 'Failed to allocate the array data'
  End If

  Do i = 1, n_data
     Write( *, * ) 'Please enter data point ', i
     Read ( *, * ) data( i )
  End Do

  arithmetic_mean =  calc_arithmetic_mean( data )
  geometric_mean  =  calc_geometric_mean ( data )
  harmonic_mean   =  calc_harmonic_mean  ( data )

  Write( *, * ) 'The arithmetic mean is: ', arithmetic_mean
  Write( *, * ) 'The geometric  mean is: ', geometric_mean
  Write( *, * ) 'The harmonic   mean is: ', harmonic_mean

  Deallocate( data )

Contains

  Real Function calc_arithmetic_mean( data )

    Real, Dimension( : ), Intent( In ) :: data

    calc_arithmetic_mean = Sum( data ) / Size( data )

  End Function calc_arithmetic_mean
  
  Real Function calc_geometric_mean( data )

    Real, Dimension( : ), Intent( In ) :: data

    calc_geometric_mean = Product( data ) ** ( 1.0 / Real( Size( data ) ) )

  End Function calc_geometric_mean

  Real Function calc_harmonic_mean( data )

    Real, Dimension( : ), Intent( In ) :: data

    calc_harmonic_mean = Real( Size( data ) ) / Sum( 1.0 / data ) 

  End Function calc_harmonic_mean

End Program means

