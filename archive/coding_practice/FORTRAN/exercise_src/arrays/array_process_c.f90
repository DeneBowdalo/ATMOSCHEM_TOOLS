PROGRAM array_process
IMPLICIT None

integer, parameter :: wp = selected_real_kind(6,35)
real (wp), dimension(9,9)  :: array, new_array
integer :: x, y

do y = 1,9
    do x = 1,9
        array(y,x) = x 
    enddo
enddo

where(mod(array/2)==1) 
    
elsewhere
    new_array=0

end where

print*, array(:,:)

END PROGRAM array_process
