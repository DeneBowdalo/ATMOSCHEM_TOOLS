PROGRAM Prime_fact
   IMPLICIT NONE
   INTEGER      :: input_integer, i
!
! Read in a number. Note the use of the do loop
! to enforce compilance with the > 2 condition
!
   DO
      WRITE(*,*)'Input an integer > 2'
      READ(*,*)input_integer
      IF(input_integer.GT.2)EXIT
      WRITE(*,*)'Error. Integer should be >2. Please try again'
   END DO
!
! remove factors of 2
!
   DO
      IF(MOD(input_integer,2).NE.0)EXIT
      input_integer = input_integer/2
      WRITE(*,*)'factor = ',2
   END DO
!
! now check odd numbers. If you have not looked at 
! Prime_a you should work out the logic of the
! next line
!
   DO i = 3, INT(SQRT(REAL(input_integer)))+1 ,2
      DO 
         IF(MOD(input_integer,i).NE.0)EXIT
         input_integer = input_integer/i
         WRITE(*,*)'factor = ',i
      END DO
   END DO

!
! if the final prime factor is greater than the square root
! of the number then it will be left at the end
!
   IF(input_integer.NE.1) WRITE(*,*) 'factor = ',input_integer
   	
END PROGRAM Prime_fact
