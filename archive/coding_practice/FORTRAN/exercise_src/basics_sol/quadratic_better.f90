PROGRAM quadratic_better
   IMPLICIT NONE 
   Integer, Parameter :: wp = Selected_real_kind( 6, 35 )
   REAL( wp ) :: a,b,c, factor, q, root1, root2
   
!
! Input a,b,c
!
   write(*,*)'Input a, b, c'
   read(*,*)a,b,c
!
! set up factor. 
!
   factor = sqrt(b*b - 4.0_wp * a * c)
!
! set up q
!
   q = - 0.5_wp * ( b + Sign( 1.0_wp, b ) * factor )

!
! compute roots
!
   root1 = q / a
   root2 = c / q
!
! Ouput result
!
   write(*,*)'The equation has roots ',root1,root2
END PROGRAM Quadratic_better
