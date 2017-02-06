PROGRAM Prime_Factors
   IMPLICIT NONE 
   integer, parameter                  :: upper_bound = 100000
   integer, parameter                  :: array_dimension = 10000
   integer                             :: n_prime_found,i,j,input_value
   integer, dimension(array_dimension) :: prime_numbers
   logical                             :: prime

!
! Store 1 and 2 as we know they are prime
!
   n_prime_found = 2
   prime_numbers(1) = 1
   prime_numbers(2) = 2
!
! Now generate all primes up to the maximum given by upper_bound
!
   do j = 3, upper_bound
      prime = .true.
      do i = 2, int(sqrt(real(j)))+1
         if(mod(j,i).eq.0)then
             prime = .false.
             exit
         end if
      end do
      if(prime)then
!
! If the number is prime increase the number found and
! store the new number. Note - we will check the array
! bounds before storing and exit with an appropriate error
! message if we move out of bounds
!
          n_prime_found = n_prime_found + 1
          if(n_prime_found.gt.array_dimension)then
             write(*,*)'Please increase the array dimensions'
             write(*,*)'I have so far found',n_prime_found,'primes'
             stop
          end if
          prime_numbers(n_prime_found) = j
      end if
   end do
!
! Input the test number
!
   do
      write(*,*)'Please input an integer >2'
      read(*,*)input_value
!
! Remember that we can only work up to numbers that are less
! than the biggest prime found squared
!
     if(input_value.gt.upper_bound*2)then
          write(*,*)'I cant work with a number that big'
          cycle
      end if
!
! Standard error check for a number > 2
!
      if(input_value.gt.2)then
          exit
      else
         write(*,*)'Error - please try again with a number > 2'
      end if
   end do
   do i = 2, n_prime_found
      if(prime_numbers(i).gt.int(real(input_value))+1)exit
      do
         if(mod(input_value,prime_numbers(i)).ne.0)then
            exit
         else
            input_value = input_value/prime_numbers(i)
            write(*,*)prime_numbers(i)
         end if
      end do
   end do
END PROGRAM Prime_Factors
