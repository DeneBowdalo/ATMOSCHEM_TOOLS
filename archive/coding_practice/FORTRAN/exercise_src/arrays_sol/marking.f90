PROGRAM Exam_Marking
   IMPLICIT NONE 
   integer, parameter                    :: number_candidates = 1000
   integer                               :: number_entered, menu_item, temp
   integer                               :: i,j,temp1
   integer, dimension(1)                 :: location
   real, dimension(3)                    :: temp_marks
   integer, dimension(number_candidates) :: candidate_number, int_average
   real, dimension(number_candidates,3)  :: candidate_marks
   REAL         :: average
   logical      :: found
   
   number_entered = 0
!
! Menu items
!
   do
      write(*,*)'Enter 1 to enter candidate data'
      write(*,*)'Enter 2 to print a candidates marks'
      write(*,*)'Enter 3 to print a candidates average and letter grade'
      write(*,*)'Enter 4 to print the class list'
      write(*,*)'Enter 5 to exit'
      read(*,*)menu_item
!
! First check - if we have entered no candidates then 2-4 all invalid
!
      if((number_entered.eq.0).and.(menu_item.ne.1))then
          write(*,*)'Please enter a candidates data before any other operation'
          cycle
      end if
!
! Second check - item must be 1-5 inclusive
!
      if(menu_item.lt.1.or.menu_item.gt.5)then
          write(*,*)'Only selections between 1 and 5 are valid'
          cycle
      end if
!
! Exit condition
!
      if(menu_item.eq.5)exit
!
! Menu items chosen via select case
! 
      select case (menu_item)
         case(1)  ! Enter data
            write(*,*)'Input candidate number.'
            read(*,*)temp
!
! Bit of logic to see if we are trying to overwrite an existing candidadate
!
            if(number_entered.gt.0)then
                found = .false.
                do i = 1, number_entered
                   if(temp.eq.candidate_number(i))then
                       found = .true.
                       temp1 = i
                       exit
                   end if
                end do
                if(found)then
                    write(*,*)'This candidate is already on the list'
                    write(*,*)'Do you wish to overwrite? 1 for yes 2 for no'
                    read(*,*)temp
                    if(temp.ne.1)cycle
                else
                   number_entered = number_entered + 1
                   temp1 = number_entered
                   candidate_number(number_entered) = temp
                end if
            else
                number_entered = 1
                temp1 = number_entered
                candidate_number(number_entered) = temp
            end if
            write(*,*)'Input three marks'
            read(*,*)(candidate_marks(temp1,i),i=1,3)
            average = sum(candidate_marks(temp1,:))/3
            int_average(temp1) = int(average)
            if(average - real(int_average(temp1)).gt.0.5) &
               int_average(temp1) = int_average(temp1) + 1
            cycle
         case(2)
             write(*,*)'Input candidate number.'
             read(*,*)temp
             found = .false.
             do i = 1, number_entered
                if(temp.eq.candidate_number(i))then
                    write(*,*)candidate_marks(i,:)
                    found = .true.
                    exit
                end if
             end do
             if(.not.found)then
                 write(*,*)'Candidate ',temp,' has not been entered'
             end if
             cycle
         case(3)
             write(*,*)'Input candidate number.'
             read(*,*)temp
             found = .false.
             do i = 1, number_entered
                if(temp.eq.candidate_number(i))then
                    found = .true.
                    temp = i
                    exit
                end if
             end do
             if(.not.found)then
                 write(*,*)'Candidate ',temp,' has not been entered'
             else
                 write(*,*)'Average Mark =',int_average(temp)
                 select case(int_average(temp))
                    case(90:)
                       write(*,*)'Grade A'
                    case(85:89)
                       write(*,*)'Grade AB'
                    case(80:84)
                       write(*,*)'Grade B'
                    case(75:79)
                       write(*,*)'Grade BC'
                    case(70:74)
                       write(*,*)'Grade C'
                    case(65:69)
                       write(*,*)'Grade CD'
                    case(60:64)
                       write(*,*)'Grade D'
                    case(:59)
                       write(*,*)'Grade F'
                 end select
             end if
             cycle
         case(4)
              do i = 1, number_entered-1
                  location = maxloc(int_average(i:number_entered))
                  j = location(1) + i - 1
                  temp = int_average(j)
                  int_average(j) = int_average(i)
                  int_average(i) = temp
                  temp = candidate_number(j)
                  candidate_number(j) = candidate_number(i)
                  candidate_number(i) = temp
                  temp_marks = candidate_marks(j,:)
                  candidate_marks(j,:) = candidate_marks(i,:)
                  candidate_marks(i,:) = temp_marks
               end do
               do i = 1, number_entered
                  write(*,*)candidate_number(i),int_average(i)
               end do
               cycle
      end select
   end do
END PROGRAM Exam_Marking
