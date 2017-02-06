Module io
   Implicit None
   Public :: write
   Private !set everything to private unless otherwise defined
!
!use an interface to create an overloaded procedure
!
   Interface write
      Module Procedure write_integer_1d
      Module Procedure write_integer_2d
   End Interface
Contains
   Subroutine write_integer_1d( i )
      Integer, Dimension( : ), Intent( In ) :: i
      Write( *, * ) i
   End Subroutine write_integer_1d
   Subroutine write_integer_2d( i )
      Integer, Dimension( :, : ), Intent( In ) :: i
      Write( *, * ) i
   End Subroutine write_integer_2d
End Module io
