PROGRAM means
IMPLICIT None

Real, Dimension( : ), Allocatable :: data
character (len=1) :: type_mean
integer :: n_num, error, i
real :: mean
do

print*, 'How many numbers would you like to calculate the mean of?'
read(*,*) n_num

Allocate( data( 1:n_num ), Stat = error )
  If( error /= 0 ) Then
     Write( *, * ) 'Failed to allocate the array data'
  End If

  Do i = 1, n_num
     Write( *, * ) 'Please enter data point ', i
     Read ( *, * ) data( i )
  End Do

print*,'What kind of mean would you like to calculate?'
print*,'Type a for Arithmetic Mean, g for Geometric Mean and h for Harmonic Mean' 
read(*,*) type_mean

if(type_mean == 'a') then
    call arith(data, mean)
    print*,''
    print*,''
else if(type_mean == 'g') then
    call geometric(data, mean)
    print*,''
    print*,''
else if(type_mean == 'h') then
    call harmonic(data, mean) 
    print*,''
    print*,''
else 
    print*, 'Please enter a valid option'

endif

DEALLOCATE(data)

enddo

Contains

  Subroutine arith( data, mean )

    Real, Dimension( : ), Intent( In    ) :: data
    Real                , Intent(   Out ) :: mean

    mean = Sum( data ) / Size( data )
    print*,'Arithmetic mean equals', mean
  End Subroutine arith

  Subroutine geometric( data, mean )

    Real, Dimension( : ), Intent( In    ) :: data
    Real                , Intent(   Out ) :: mean

    mean = Product( data ) ** ( 1.0 / Real( Size( data ) ) )
    print*,'Geometric mean equals', mean
  End Subroutine geometric

  Subroutine harmonic( data, mean )

    Real, Dimension( : ), Intent( In    ) :: data
    Real                , Intent(   Out ) :: mean

    mean = Real( Size( data ) ) / Sum( 1.0 / data )
    print*,'Harmonic Mean Equals:', mean
  End Subroutine harmonic

END PROGRAM means
