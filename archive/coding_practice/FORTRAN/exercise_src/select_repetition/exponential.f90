PROGRAM Exponential 
   IMPLICIT NONE 
   INTEGER :: Counter
   REAL :: Begin, r_End, Step, X, Tolerance, ExpX, Term, Sums
   WRITE(*,*) 'input Initial, Final, Step and Tolerance' 
   READ(*,*) Begin, r_End, Step, Tolerance 
   WRITE(*,*) '     X      |     Exp(X)    |  Answer  | Number of Iterations'
   X = Begin
   DO ! Loop over x values
      IF (X > r_End) EXIT
!
!use the fortran intrinsic to get a solution to compare to our answer 
!
      ExpX = EXP(X)
      Counter = 1 
      Term = X 
      Sums = 1.0
      DO ! Do Loop until we meet required tolerence
        IF (ABS(Term) < Tolerance) EXIT
        Sums = Sums + Term
        Counter = Counter + 1
        Term = Term * (X / Counter)
      END DO
      WRITE(*,*) X,' ',ExpX,' ',Sums,' ',Counter
      X = X + Step 
   END DO 
END PROGRAM Exponential 

