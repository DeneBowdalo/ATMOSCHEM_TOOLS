PROGRAM Multiply
  IMPLICIT NONE

  INTEGER :: i, j
  INTEGER,PARAMETER :: n = 10
  REAL, DIMENSION(n,n) :: a, b, c1, c2, c3
  
  CALL RANDOM_NUMBER( a )
  CALL RANDOM_NUMBER( b )
!
! Whole array operations are component wise 
!  
 c1 = a * b

  WRITE(*,*) 'All greater than zero? ', ALL( c1>0 )
  WRITE(*,*) 'Amount greater or equal to 0.5: ', COUNT( c1>=0.5 )
  WRITE(*,*) 'Smallest value: ', MINVAL( c1 )
  WRITE(*,*) 'Largestest value: ', MAXVAL( c1 )

  WHERE ( c1 < 0.5 )
     c1 = 1.0
  END WHERE

  WRITE(*,*) 'Amount greater or equal to 0.5 now: ', COUNT( c1>=0.5 ), &
       '( Should be: ', n**2, ')'    
  WRITE(*,*)
!
! Matrix multiplication
!
  DO j = 1, n
     DO i = 1, n

        c1(i,j) = SUM( a(i,:) * b(:,j) )
        
     END DO
  END DO

!
! Do the multiplication with an intrinsic
!
  c3 = MATMUL(a,b)
!
! Check c1 and c2 are the same as c3
!
  WRITE(*,*) 'Error in do loop mult: ', SUM(c1-c3)
  WRITE(*,*)

END PROGRAM Multiply
