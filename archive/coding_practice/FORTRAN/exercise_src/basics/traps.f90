PROGRAM      Fortran_Traps
   IMPLICIT     NONE

   INTEGER, PARAMETER          :: A = 2, B = 2, H = 3
   INTEGER, PARAMETER          :: O = 4, P = 6
   CHARACTER(LEN=5), PARAMETER :: M = 'Smith', N = 'TEXAS'
   CHARACTER(LEN=4), PARAMETER :: X = 'Smith'
   CHARACTER(LEN=6), PARAMETER :: Y = 'TEXAS'

!  The exponential trap

   WRITE(*,*)   "First, the exponential trap:"
   WRITE(*,*)   A, ' ** ', B, ' ** ', H, ' = ', A**B**H
   WRITE(*,*)   '( ', A, ' ** ', B, ' ) **', H, ' = ', (A**B)**H
   WRITE(*,*)   A, ' ** ( ', B, ' ** ', H, ' ) = ', A**(B**H)
   WRITE(*,*)

!  The integer division trap.  Intrinsic function REAL() converts
!  an integer to a real number

   WRITE(*,*)   "Second, the integer division trap:"
   WRITE(*,*)
   WRITE(*,*)   O, ' / ', P, ' = ', O/P
   WRITE(*,*)   'REAL( ', O, ' ) / ', P, ' = ', REAL(O)/P
   WRITE(*,*)   O, ' / REAL( ', P, ' ) = ', O/REAL(P)
   WRITE(*,*)

!  The string truncation trap

   WRITE(*,*)   "Third, the string truncation trap:"
   WRITE(*,*)   'IS ', M, ' STILL IN ', N, '?'
   WRITE(*,*)   'IS ', X, ' STILL IN ', Y, '?'

END PROGRAM  Fortran_Traps
