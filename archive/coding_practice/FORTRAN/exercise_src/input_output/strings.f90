Program Strings
   implicit none
   integer                         :: i
   character(len=17), dimension(5) :: parts

   parts(1)='The quick brown  '
   parts(2)='fox '
   parts(3)='jumped over the '
   parts(4)='lazy'
   parts(5)='dog.'
!
! Write with single write statement
!
   write(*,*)parts
!
! Now write with do loop
!
   do i = 1, 5
      write(*,*)parts(i)
   end do
!
! Insert formatting into the next write statement
! to produce a single line of correctly formatted output
!
   write(*,*)parts
!
! Insert formatting into the next write statement to
! produce a column of output in which each line is left justified
! to the start of the words
!
   do i = 1, 5
      write(*,*)parts(i)
   end do
End Program Strings
