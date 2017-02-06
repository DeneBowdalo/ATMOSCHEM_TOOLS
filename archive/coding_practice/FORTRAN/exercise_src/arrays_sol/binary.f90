PROGRAM Binary_Search
   IMPLICIT NONE
   integer                            :: n, alloc_err, i, j, temp
   integer                            :: search_value, left, right, mid
   integer, dimension(1)              :: location
   integer, dimension(:), allocatable :: x
!
! input the number of elements in the array
!
   do
      write(*,*)'Input number of elements in array'
      read(*,*)n
      if(n.le.1)then
          write(*,*)'The array should have more than 1 element'
      else
          exit
      end if
   end do
!
! Allocate memory for the array and input
!
   allocate(x(n),stat=alloc_err)
   if(alloc_err.ne.0)then
       write(*,*)'Error allocating array'
       stop
   end if
   do i = 1, n
      write(*,*)'Input element ' ,i,' of array'
      read(*,*)x(i)
   end do
!
! Now sort
!
   do i = 1, n-1
      location = minloc(x(i:n))
!
! Note two things; firstly minloc returns an array.
! Secondly it returns the location of the minimum within the
! array segment passed to it, not within the whole array. Hence
! we need to add (i-1) to its value to take account of this offset
!
      j = location(1) + i - 1
      temp = x(j)
      x(j) = x(i)
      x(i) = temp
   end do
!
! and output the sorted array
!
   write(*,*)'The sorted array is;'
   write(*,*)x
!
! Set up left and right as the extrema
!
   left = 1
   right = n
!
! Input search number
!
   write(*,*)'Input Search Number'
   read(*,*)search_value
!
! Is search_value within the limis of the array?
!
   if(search_value.lt.x(left))then
      write(*,*)'Value is outside the left bound of the array'
   else if(search_value.gt.x(right))then
      write(*,*)'Value is outside the right bound of the array'
   else if(search_value.eq.x(left))then
      write(*,*)'Value is located at position',left
   else if(search_value.eq.x(right))then
      write(*,*)'Value is located at position',right
   else   ! value is inside the array
      do
         mid = (right + left)/2
         if(x(mid).eq.search_value)then
             write(*,*)'The value is located at position',mid
             exit
         else
             if(x(mid).lt.search_value)then
                 left = mid
             else
                 right = mid
             end if
         end if
!
! We have to be careful as left and right approach each other
!
         if(right-left.eq.1)then
             write(*,*)'The value is bracketed by positions',right,left
             exit
         end if
      end do
   end if
!
! and tidy up by deallocating memory
!
   deallocate(x)
END PROGRAM Binary_Search

