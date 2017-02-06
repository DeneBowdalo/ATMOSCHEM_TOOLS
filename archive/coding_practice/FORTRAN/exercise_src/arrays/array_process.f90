PROGRAM array_process
IMPLICIT None

integer, parameter :: wp = selected_real_kind(6,35)
real (wp), dimension(9,9)  :: array
integer :: x, y

do y = 1,9
    do x = 1,9
        array(y,x) = x 
    enddo
enddo

print*, array(1,:)

END PROGRAM array_process
