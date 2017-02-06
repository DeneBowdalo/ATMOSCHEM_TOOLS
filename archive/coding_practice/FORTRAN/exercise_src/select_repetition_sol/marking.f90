PROGRAM Marking
   IMPLICIT NONE 
   INTEGER      :: candidate_number, int_average
   REAL         :: a,b,c, average
!
! We want to repeatadly read a candidate number
! Construct a general do loop
!
   DO
!
! Read candidate number and check for exit condition
!
      WRITE(*,*)'Input candidate number. -ve number to exit'
      READ(*,*)candidate_number
      IF(candidate_number.LT.0)EXIT
!
! Read in three grades
!
      WRITE(*,*)'Input three grades'
      READ(*,*)a,b,c
!
! Compute average
!
      average = (a+b+c)/3.0
!
! And round. Note that int returns the integer part not the
! nearest integer. The formulae below will round properly. Try
! to modify the program to use the NINT intrinsic which
! returns the nearest integer rather than the integer part
!
      int_average = INT(average)
      IF(average - REAL(int_average).GT.0.5)int_average = int_average + 1
!
! We have a whole bunch of conditions which looks like a select case
! Remember that the : constructor includes = so be careful with your limits
!
      SELECT CASE(int_average)
         CASE(90:)
             WRITE(*,*)'Grade A'
         CASE(85:89)
             WRITE(*,*)'Grade AB'
         CASE(80:84)
             WRITE(*,*)'Grade B'
         CASE(75:79)
             WRITE(*,*)'Grade BC'
         CASE(70:74)
             WRITE(*,*)'Grade C'
         CASE(65:69)
             WRITE(*,*)'Grade CD'
         CASE(60:64)
             WRITE(*,*)'Grade D'
         CASE(:59)
             WRITE(*,*)'Grade F'
      END SELECT
   END DO
 END PROGRAM Marking
