PROGRAM Sorting
   IMPLICIT NONE 
   integer                         :: n, alloc_err, i, j, temp
   integer, dimension(1)           :: location
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
! Output the unsorted array
!
   write(*,*)'The unsorted array is;'
   write(*,*)x
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
! and tidy up by deallocating memory
!
   deallocate(x)
END PROGRAM Sorting
