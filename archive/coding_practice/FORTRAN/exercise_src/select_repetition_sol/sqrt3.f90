PROGRAM root
  IMPLICIT NONE
  Integer, Parameter :: wp = Selected_real_kind( 6, 35 )
  REAL( wp ) :: x0, x1
  INTEGER    :: numit, i

  PRINT *, 'Enter number of iterations'
  READ(*,*) numit
  PRINT *, 'Enter initial guess'
  READ(*,*) x0
  DO i=1,numit
     x1=.5_wp*(3.0_wp/x0+x0)
     PRINT *, x0
     x0 = x1
  END DO

END PROGRAM root
