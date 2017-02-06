PROGRAM prime

IMPLICIT NONE

integer :: try,a=2, i
logical :: primer=.TRUE.

print *,"Enter a positive number greater than two to check if it is a prime"
read(*,*) try


do i =a, try
    if(mod(try,a) == 0) then
    primer = .FALSE.
    end if
end do


if(primer == .FALSE.) then
    print *, try, "is not a prime."

else if(primer == .TRUE.) then
    print *, try, "is a prime." 

end if

END PROGRAM prime 
