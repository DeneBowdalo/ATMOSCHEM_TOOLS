PROGRAM prime
IMPLICIT NONE

integer :: try,a=3,j, i, primes
logical :: primer=.TRUE.

print *,"Enter a positive number greater than two to check all the prime numbers in the range to the number entered."
read(*,*) try

write(*,*) "The Primes from 2 to", try,"are:"


do i = a, try-1    
primer =.TRUE.
    do j = a, i-1
        primes=mod(i,j)
        if(primes == 0) then
            primer =.FALSE.
        end if
    end do
if(primer ==.TRUE.) then
    write(*,*), i
end if
end do

END PROGRAM prime
