PROGRAM numbers_test
   IMPLICIT NONE 

   Integer, Parameter :: dp = Selected_real_kind( 12, 70 )

   character        :: a

   integer          :: i = 789

   logical          :: l

   complex          :: c
   real             :: x = 123.56

   Real( dp )       :: d = 123.56d0
   Complex( dp )    :: dc
   
   write(*,*)'Kinds'
   write(*,*)'*****'
   write(*,*)'Character kind =',kind(a)
   write(*,*)'Complex kind =',kind(c)
   write(*,*)'Double precision kind=',kind(d)
   write(*,*)'Integer kind =',kind(i)
   write(*,*)'Logical kind =',kind(l)
   write(*,*)'Real kind =',kind(x)
   write(*,*)
   write(*,*)'Integers:'
   write(*,*)'*********'
   !
   ! radix returns the base of the model used to store integers.
   ! Usually we will be working in binary and radix should return 2.
   ! bit_size returns the total number of bits used to store integers.
   ! digits returns the total number of bits used to represent the
   ! magnitude of the integer. Usually this is bit_size - 1. The extra
   ! bit is the sign bit. The largest number should be radix**digits -1
   !
   write(*,*)'Base of model =',radix(i)
   write(*,*)'Bit Size =',bit_size(i)
   write(*,*)'Significant bits =',digits(i)
   write(*,*)'Max value of i =',huge(i)
   write(*,*)
   write(*,*)'Reals:'
   write(*,*)'******'
   write(*,*)'Base of model',radix(x)
   write(*,*)'Max value of x =',huge(x)
   write(*,*)'Min value of x =',tiny(x)
   write(*,*)'Decimal precision =',precision(x)
   write(*,*)'Negligible value =',epsilon(x)
   write(*,*)'Nearest +ve =',nearest(x,1.0)
   write(*,*)'Nearest -ve =',nearest(x,-1.0)
   write(*,*)'Absolute spacing =',spacing(x)
   write(*,*)
   write(*,*)'"Double":'
   write(*,*)'*****************'
   write(*,*)'Base of model',radix(d)
   write(*,*)'Max value of x =',huge(d)
   write(*,*)'Min value of x =',tiny(d)
   write(*,*)'Decimal precision =',precision(d)
   write(*,*)'Negligible value =',epsilon(d)
   write(*,*)'Nearest +ve =',nearest(d,1.0)
   write(*,*)'Nearest -ve =',nearest(d,-1.0)
   write(*,*)'Absolute spacing =',spacing(d)
   write(*,*)
   write(*,*)'Complex'
   write(*,*)'*******'
   write(*,*)'Decimal precision =',precision(c)
   write(*,*)
   write(*,*)'"Double" Complex'
   write(*,*)'*******'
   write(*,*)'Decimal precision =',precision(dc)
   write(*,*)
END PROGRAM numbers_test
