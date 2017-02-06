PROGRAM Array_c
   IMPLICIT NONE
   integer                 :: i,j
   integer, dimension(9,9) :: array
!
!loop over rows and columns to assign elements of the array
!remember Fortran has column major order so loop over j then i 
!   
   do i = 1,9
      do j = 1,9
         array(j,i) = i+10*j
      end do
   end do
!
!set all the even values (those that have a remainder of 0 when
!divided by 2) to 0
!
   where(mod(array,2).eq.0)
       array = 0
   end where
   write(*,*)array
END PROGRAM Array_c
