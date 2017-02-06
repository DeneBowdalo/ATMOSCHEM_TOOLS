PROGRAM quadratic
   IMPLICIT NONE 
   Integer, Parameter :: wp = Selected_real_kind( 12, 70 )
   REAL( wp ) :: a,b,c, factor, root1, root2
   
!
! Input a,b,c
!
   write(*,*)'Input a, b, c'
   read(*,*)a,b,c
!
! set up factor. Used in root1 and root2 so compute once
!
   factor = sqrt(b*b - 4.0_wp * a * c)
!
! compute roots
!
   root1 = (-b + factor)/(2.0_wp*a)
   root2 = (-b - factor)/(2.0_wp*a)
!
! Ouput result
!
   write(*,*)'The equation has roots ',root1,root2
END PROGRAM quadratic

