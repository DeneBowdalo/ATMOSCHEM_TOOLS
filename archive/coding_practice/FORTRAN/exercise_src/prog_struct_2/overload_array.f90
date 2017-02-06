Program overload_array
   Use io
   Implicit None
!
!declare and initialise two arrays of different rank
!
   Integer, Dimension(      1:8 ) :: i = -5
   Integer, Dimension( 1:3, 1:3 ) :: j = 3
!
!call the overloaded subroutine write to print the arrays
!
   Call write( i )
   Call write( j )
End Program overload_array
