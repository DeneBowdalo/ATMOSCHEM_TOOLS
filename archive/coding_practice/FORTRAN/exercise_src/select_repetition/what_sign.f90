PROGRAM what_sign

Implicit None
integer :: intn

do while (intn/=99999)  

    read(*,*) intn

    if(intn > 0) then
        write (*,*) intn, 'is a positive number'

    else if(intn < 0) then
        write (*,*) intn, 'is a negative number'

    else
        write (*,*) intn, 'is equal to zero'

    end if

end do
END PROGRAM
