Program Strings
   IMPLICIT NONE
   integer                         :: i
   character(len=17), dimension(5) :: parts
   character(len=8), dimension(5) :: format_strings

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
   write(*,'(a16,a4,a16,a4,1x,a4)')parts
!
! Insert formatting into the next write statement to
! produce a column of output in which each line is left justified
! to the start of the words
!
   format_strings(1) = '(a15)'
   format_strings(2) = '(12x,a3)'
   format_strings(3) = format_strings(1)
   format_strings(4) = '(11x,a4)'
   format_strings(5) = format_strings(4)
   do i = 1, 5
      write(*,format_strings(i))parts(i)
   end do
End Program Strings
