PROGRAM means
   IMPLICIT NONE 
   REAL :: a,b,c
!
! input a,b and c
!
   write(*,*)'Please input three numbers'
   read(*,*)a,b,c
!
! compute means and output
!
   write(*,*)'Arithmetic mean =',(a+b+c)/3.0
   write(*,*)'Geometric mean =',(a*b*c)**(0.333333333)
   write(*,*)'Harmonic mean =',3.0/(1.0/a + 1.0/b + 1.0/c)
END PROGRAM Means
