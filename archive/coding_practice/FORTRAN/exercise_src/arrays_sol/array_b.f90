PROGRAM Array_b
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
!print the whole of the 1st row
!
   write(*,*)array(:,1)
END PROGRAM Array_b
